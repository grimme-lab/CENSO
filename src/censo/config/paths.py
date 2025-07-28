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

from ..params import PLENGTH, DIGILEN


class PathsConfig(BaseModel):
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
    def __str__(self):
        """
        Print out the paths of all external QM programs.
        """
        # Create an empty list to store the lines of the output.
        lines = []

        # Append a separator line to the output.
        lines.append("\n" + "".ljust(PLENGTH, "-") + "\n")

        # Append the title of the section to the output, centered.
        lines.append("PATHS of external QM programs".center(PLENGTH, " ") + "\n")

        # Append a separator line to the output.
        lines.append("".ljust(PLENGTH, "-") + "\n")

        # Iterate over each program and its path in the settings.
        for program, path in self:
            # Append a line with the program and its path to the output.
            lines.append(f"{program}:".ljust(DIGILEN, " ") + f"{path}\n")

        return "".join(lines)

    @property
    def orcaversion(self) -> str:
        return self._orcaversion

    @field_validator("orca")
    def validate_orca(cls, value: str):
        if not Path(value).is_file():
            raise ValueError(f"orca executable not found at {value}.")
        return value

    @field_validator("tm")
    def validate_turbomole(cls, value: str):
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
        if not Path(value).is_file():
            raise ValueError(f"xtb executable not found at {value}.")
        return value

    @field_validator("cosmotherm")
    def validate_cosmotherm(cls, value: str):
        if not Path(value).is_file():
            raise ValueError(f"cosmotherm executable not found at {value}.")
        return value

    @field_validator("cosmorssetup")
    def validate_cosmors(cls, value: str, info: ValidationInfo):
        """
        Validate cosmotherm and cosmorssetup.
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
        # if orca was found try to determine orca version from the path (kinda hacky)
        if self.orca:
            match = re.search(r"(\d+\.\d+\.\d+)", str(self.orca))
            if match:
                self._orcaversion = match.group(1)
            else:
                # Try to extract version from binary content
                with open(self.orca, "rb") as f:
                    binary_content = f.read()

                version_pattern = rb"Program Version (\d+\.\d+\.\d+)"
                match = re.search(version_pattern, binary_content)

                if not match:
                    raise ValueError("Could not determine ORCA version.")
                else:
                    version_bytes = match.group(1)
                    version_string = version_bytes.decode("utf-8")
                    self._orcaversion = version_string

        return self
