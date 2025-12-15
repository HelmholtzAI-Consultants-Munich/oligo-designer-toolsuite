from os import PathLike
from effidict import EffiDict, LRUReplacement, PickleBackend

from .interface import KeyStore


class EffiDictKeyStore(KeyStore):
    def __init__(self, storage_path: PathLike, max_entries_in_memory: int):
        backend = PickleBackend(storage_path=storage_path)
        strategy = LRUReplacement(
            disk_backend=backend,
            max_in_memory=max_entries_in_memory,
        )
        self._store = EffiDict(
            disk_backend=backend,
            replacement_strategy=strategy,
        )

    def __getitem__(self, key: str):
        return self._store[key]

    def __setitem__(self, key: str, value):
        self._store[key] = value

    def __delitem__(self, key: str):
        del self._store[key]

    def __contains__(self, key: str) -> bool:
        return key in self._store

    def keys(self):
        return self._store.keys()

    def items(self):
        return self._store.items()
