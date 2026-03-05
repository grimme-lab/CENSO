from typing import override
import re
from pydantic import (
    BaseModel,
    Field,
    PrivateAttr,
    field_validator,
    ConfigDict,
    ValidationInfo,
    model_validator,
)
from pathlib import Path
import mmap

from ..params import PLENGTH
from ..utilities import h2


class PathsConfig(BaseModel):
    """
    Configuration for paths to external programs.
    """

    model_config = ConfigDict(
        str_strip_whitespace=True,
        use_attribute_docstrings=True,
        validate_assignment=True,
    )

    orca: str = Field("")
    """Absolute path to the orca binary."""

    _orcaversion: str = PrivateAttr("")
    """ORCA version string. Should be extracted from somewhere."""

    tm: str = Field("")
    """Absolute path to the turbomole binary directory."""

    xtb: str = Field("")
    """Absolute path to the xtb binary."""

    cosmotherm: str = Field("")
    """Absolute path to the cosmotherm binary."""

    cosmorssetup: str = Field("")
    """Name of the cosmors parameterization file."""

    @override
    def __str__(self) -> str:
        lines: list[str] = []
        lines.append(h2("PATHS of External Programs"))
        kv = [(k, v) for k, v in self]
        kv.sort(key=lambda x: type(x[1]).__name__)
        for name, value in kv:
            lines.append(f"{name:>{PLENGTH // 2 - 2}} : {value}")

        return str("\n".join(lines))

    @property
    def orcaversion(self) -> str:
        """
        Get the ORCA version string.

        :return: The ORCA version.
        """
        return self._orcaversion

    @field_validator("orca")
    def validate_orca(cls, value: str):
        """
        Validate the ORCA executable path.

        :param value: The path to validate.
        :return: The validated path.
        """
        if not Path(value).is_file():
            raise ValueError(f"orca executable not found at {value}.")
        return value

    @field_validator("tm")
    def validate_turbomole(cls, value: str):
        """
        Validate the Turbomole binary directory.

        :param value: The path to validate.
        :return: The validated path.
        """
        tmp = Path(value)
        if not all(
            (
                (tmp / "ridft").is_file(),
                (tmp / "mpshift").is_file(),
                (tmp / "escf").is_file(),
            )
        ):
            raise ValueError(
                f"Turbomole path does not contain either ridft, mpshift or escf at {value}."
            )
        return value

    @field_validator("xtb")
    def validate_xtb(cls, value: str):
        """
        Validate the XTB executable path.

        :param value: The path to validate.
        :return: The validated path.
        """
        if not Path(value).is_file():
            raise ValueError(f"xtb executable not found at {value}.")
        return value

    @field_validator("cosmotherm")
    def validate_cosmotherm(cls, value: str):
        """
        Validate the COSMOtherm executable path.

        :param value: The path to validate.
        :return: The validated path.
        """
        if not Path(value).is_file():
            raise ValueError(f"cosmotherm executable not found at {value}.")
        return value

    @field_validator("cosmorssetup")
    def validate_cosmors(cls, value: str, info: ValidationInfo):
        """
        Validate cosmotherm and cosmorssetup.

        :param value: The cosmorssetup value.
        :param info: Validation info.
        :return: The validated value.
        """
        if "cosmotherm" in info.data and info.data["cosmotherm"]:
            cosmotherm_path = Path(info.data["cosmotherm"])
            # Re-validate cosmotherm in case it was modified
            if not cosmotherm_path.is_file():
                raise ValueError(
                    f"cosmotherm executable not found at {cosmotherm_path}."
                )

            # Now validate cosmorssetup
            if value:
                setupfile = (
                    cosmotherm_path.parent / ".." / "CTDATA-FILES" / value
                ).resolve()
                if not setupfile.is_file():
                    raise ValueError(f"{setupfile} not found.")
        return value

    @model_validator(mode="after")
    def get_orca_version(self):
        """
        Extract and set the ORCA version.

        :return: The validated instance.
        """
        if self.orca:
            # Try parsing version from binary directly
            with open(self.orca, "rb") as f:
                version_pattern = rb"Program Version (\d+\.\d+\.\d+)"
                version_bytes: bytes | None = None

                try:
                    with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                        match_bytes = re.search(version_pattern, mm)
                        if match_bytes:
                            # Copy the matched bytes out of the mmap buffer before
                            # it is closed; accessing group(1) after context exit
                            # raises TypeError on some platforms.
                            version_bytes = bytes(match_bytes.group(1))
                except (OSError, ValueError):
                    # Some files (or mocked files in tests) cannot be memory-mapped,
                    # e.g. empty files. Fall back to reading bytes directly.
                    f.seek(0)
                    match_bytes = re.search(version_pattern, f.read())
                    if match_bytes:
                        version_bytes = bytes(match_bytes.group(1))

                if not version_bytes:
                    raise ValueError(
                        f"Could not determine ORCA version. Please check {self.orca}"
                    )

                self._orcaversion = version_bytes.decode("utf-8")

        return self
