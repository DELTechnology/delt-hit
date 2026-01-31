from pydantic import BaseModel

class Region(BaseModel):
    name: str
    index: int
    codons: list[str]
    max_error_rate: float
    indels: int

    @property
    def id(self):
        return f'{self.index}-{self.name}'
