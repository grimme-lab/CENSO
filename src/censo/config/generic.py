import enum
from pydantic import BaseModel, ConfigDict, field_serializer


class GenericConfig(BaseModel):
    """
    Generic configuration base class using Pydantic for validation and serialization.
    """

    model_config = ConfigDict(
        # use_enum_values=True,  # this might lead to problems somewhere
        str_to_lower=True,
        str_strip_whitespace=True,
        validate_default=True,
        use_attribute_docstrings=True,  # , validate_assignment=True <- don't use this for now
    )

    @field_serializer("*", when_used="json-unless-none")
    def serialize_enums(self, value):
        """
        Serialize enum values to their string representation for JSON output.

        :param value: The value to serialize.
        :return: The serialized value, with enums converted to their values.
        """
        if isinstance(value, enum.Enum):
            return value.value
        return value
