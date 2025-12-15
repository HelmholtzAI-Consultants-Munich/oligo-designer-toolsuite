
from os import PathLike
from typing import Literal, TypeAlias

from .interface import KeyStore


KeystoreKind = Literal['effidict', 'inmemory', 'msgspec', 'pickle', 'sqlite']


def make_keystore(kind: KeystoreKind, path: PathLike, inmemory_cache_size: int = 0) -> KeyStore:
    if kind == 'effidict':
        from .effidict_keystore import EffiDictKeyStore
        return EffiDictKeyStore(storage_path=path, max_entries_in_memory=inmemory_cache_size)
    elif kind == 'inmemory':
        from .inmemory_keystore import InMemoryKeyStore
        return InMemoryKeyStore()
    elif kind == 'msgspec':
        from .msgspec_keystore import MsgspecKeyStore
        return MsgspecKeyStore(base_dir=path)
    elif kind == 'pickle':
        from .pickle_keystore import PickleKeyStore
        return PickleKeyStore(base_dir=path)
    elif kind=='sqlite':
        from .sqlite_keystore import SQLiteKeyStore
        return SQLiteKeyStore(db_path=path)
    
    else:
        raise ValueError(f'keystore {kind} not found.')