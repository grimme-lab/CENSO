import enum
from pydantic import BaseModel, ConfigDict, field_serializer


class GenericConfig(BaseModel):
    model_config = ConfigDict(
        # use_enum_values=True,  # this might lead to problems somewhere
        str_to_lower=True,
        str_strip_whitespace=True,
        validate_default=True,
        use_attribute_docstrings=True,  # , validate_assignment=True <- don't use this for now
    )

    @field_serializer("*", when_used="json-unless-none")
    def serialize_enums(self, value):
        """If the value is an enum, return its value, otherwise return as is."""
        if isinstance(value, enum.Enum):
            return value.value
        return value
