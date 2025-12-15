from os import PathLike
from .interface import KeyStore


class InMemoryKeyStore(KeyStore):
    def __init__(self):
        self._store: dict[str, dict] = {}

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
