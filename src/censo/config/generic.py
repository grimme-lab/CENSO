from pydantic import BaseModel, ConfigDict


class GenericConfig(BaseModel):
    model_config = ConfigDict(
        use_enum_values=True,  # this might lead to problems somewhere
        str_to_lower=True,
        str_strip_whitespace=True,
        validate_default=True,
        use_attribute_docstrings=True,  # , validate_assignment=True <- don't use this for now
    )
