from pydantic import BaseModel

class Region(BaseModel):
    """Validated region definition for demultiplexing."""
    name: str
    index: int
    codons: list[str]
    max_error_rate: float
    indels: int

    @property
    def id(self):
        """Return a stable region identifier."""
        return f'{self.index}-{self.name}'
