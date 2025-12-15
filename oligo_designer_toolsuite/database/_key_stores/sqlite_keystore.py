from os import PathLike
import sqlite3
from typing import Any
import msgspec
import pandas as pd
from pathlib import Path

from .interface import KeyStore


class SQLiteKeyStore(KeyStore):
    """
    SQLite-backed region store.
    Regions stored as msgspec-encoded blobs.
    """

    def __init__(self, db_path: PathLike):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)

        self._encoder = msgspec.msgpack.Encoder()
        self._decoder = msgspec.msgpack.Decoder()

        self._conn = sqlite3.connect(self.db_path)
        self._init_db()

    def _init_db(self):
        cur = self._conn.cursor()
        cur.execute(
            """
            CREATE TABLE IF NOT EXISTS regions (
                region_id TEXT PRIMARY KEY,
                payload   BLOB NOT NULL
            )
            """
        )
        self._conn.commit()

    def __contains__(self, key: str) -> bool:
        cur = self._conn.cursor()
        cur.execute("SELECT 1 FROM regions WHERE region_id = ?", (key,))
        return cur.fetchone() is not None

    def __getitem__(self, key: str) -> Any:
        cur = self._conn.cursor()
        cur.execute(
            "SELECT payload FROM regions WHERE region_id = ?", (key,)
        )
        row = cur.fetchone()
        if row is None:
            raise KeyError(key)

        data = self._decoder.decode(row[0])
        return self._post_load(data)

    def __setitem__(self, key: str, value: Any) -> None:
        data = self._encoder.encode(self._pre_dump(value))
        cur = self._conn.cursor()
        cur.execute(
            """
            INSERT INTO regions (region_id, payload)
            VALUES (?, ?)
            ON CONFLICT(region_id)
            DO UPDATE SET payload = excluded.payload
            """,
            (key, data),
        )
        self._conn.commit()

    def __delitem__(self, key: str) -> None:
        cur = self._conn.cursor()
        cur.execute(
            "DELETE FROM regions WHERE region_id = ?", (key,)
        )
        self._conn.commit()

    def keys(self):
        cur = self._conn.cursor()
        cur.execute("SELECT region_id FROM regions")
        return (row[0] for row in cur.fetchall())

    def items(self):
        for key in self.keys():
            yield key, self[key]

    # ---- helpers ----

    def _pre_dump(self, obj: Any) -> Any:
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

    def close(self):
        self._conn.close()
