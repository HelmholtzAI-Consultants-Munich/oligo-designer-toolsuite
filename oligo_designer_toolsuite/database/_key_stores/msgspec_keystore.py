from os import PathLike
from pathlib import Path
from typing import Any
import msgspec
import pandas as pd

from .interface import KeyStore

class MsgspecKeyStore(KeyStore):
    """
    Explicit per-region storage using msgspec + msgpack.
    One file per region.
    """

    def __init__(self, base_dir: PathLike):
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)

        self._encoder = msgspec.msgpack.Encoder()
        self._decoder = msgspec.msgpack.Decoder()

    def _path(self, region_id: str) -> Path:
        return self.base_dir / f"{region_id}.msgpack"

    def __contains__(self, key: str) -> bool:
        return self._path(key).exists()

    def __getitem__(self, key: str) -> Any:
        path = self._path(key)
        if not path.exists():
            raise KeyError(key)

        with open(path, "rb") as f:
            data = self._decoder.decode(f.read())

        return self._post_load(data)

    def __setitem__(self, key: str, value: Any) -> None:
        path = self._path(key)
        data = self._pre_dump(value)

        with open(path, "wb") as f:
            f.write(self._encoder.encode(data))

    def __delitem__(self, key: str) -> None:
        path = self._path(key)
        if path.exists():
            path.unlink()

    def keys(self):
        return (p.stem for p in self.base_dir.glob("*.msgpack"))

    def items(self):
        for key in self.keys():
            yield key, self[key]

    # ---- helpers ----

    def _pre_dump(self, obj: Any) -> Any:
        """Convert non-msgpack-safe objects."""
        if isinstance(obj, pd.DataFrame):
            return {
                "__type__": "pandas.DataFrame",
                "value": obj.to_dict(orient="split"),
            }
        return obj

    def _post_load(self, obj: Any) -> Any:
        if isinstance(obj, dict) and obj.get("__type__") == "pandas.DataFrame":
            return pd.DataFrame(**obj["value"])
        return obj
