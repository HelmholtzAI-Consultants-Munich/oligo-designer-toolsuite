from os import PathLike
from pathlib import Path
from typing import Any
import pickle

from .interface import KeyStore


class PickleKeyStore:
    """
    Explicit per-region storage using pickle.
    One file per region.
    """

    def __init__(self, base_dir: PathLike):
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)

    def _path(self, key: str) -> Path:
        return self.base_dir / f"{key}.pkl"

    def __contains__(self, key: str) -> bool:
        return self._path(key).exists()

    def __getitem__(self, key: str) -> Any:
        path = self._path(key)
        if not path.exists():
            raise KeyError(key)

        with open(path, "rb") as f:
            return pickle.load(f)

    def __setitem__(self, key: str, value: Any) -> None:
        path = self._path(key)
        with open(path, "wb") as f:
            pickle.dump(value, f, protocol=pickle.HIGHEST_PROTOCOL)

    def __delitem__(self, key: str) -> None:
        path = self._path(key)
        if path.exists():
            path.unlink()

    def keys(self):
        return (p.stem for p in self.base_dir.glob("*.pkl"))

    def items(self):
        for key in self.keys():
            yield key, self[key]
