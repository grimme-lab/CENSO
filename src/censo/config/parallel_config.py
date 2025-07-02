from pydantic import BaseModel, Field, model_validator

from ..params import OMPMIN_DEFAULT, OMPMAX_DEFAULT


class ParallelConfig(BaseModel):
    ncores: int = Field(gt=0)
    omp: int = Field(gt=0)
    ompmin: int = Field(gt=0, default=OMPMIN_DEFAULT)
    ompmax: int = Field(gt=0, default=OMPMAX_DEFAULT)

    @model_validator(mode="after")
    def ompmax_ge_ompmax(self):
        if self.ompmax < self.ompmin:
            raise ValueError("ompmax has to be greater or equal to ompmin.")
        return self
