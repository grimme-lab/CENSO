import json
from pathlib import Path

from ..params import USER_ASSETS_PATH

func_path = Path(__file__).parent / "dfa.json"
user_func_path = Path(USER_ASSETS_PATH) / "dfa.json"

solv_path = Path(__file__).parent / "solvents.json"
user_solv_path = Path(USER_ASSETS_PATH) / "solvents.json"

FUNCTIONALS: dict[str, dict[str, str]] = json.loads(func_path.read_text()) | (
    json.loads(user_func_path.read_text()) if user_func_path.is_file() else {}
)
SOLVENTS: dict[str, dict[str, str]] = json.loads(solv_path.read_text()) | (
    json.loads(user_solv_path.read_text()) if user_solv_path.is_file() else {}
)
