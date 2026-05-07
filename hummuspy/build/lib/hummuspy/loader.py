from importlib.resources import files

tfs_resources = [
    "human_tfs_r_hummus",
    "mouse_tfs_r_hummus",
]


def load_tfs(tfs_list) -> list[str]:
    if tfs_list not in tfs_resources:
       raise ValueError(f"The name of the tfs list must be in {tfs_resources}")

    text = files("hummuspy.data").joinpath(tfs_list+".txt").read_text(encoding="utf-8")
    return [line for line in (s.strip() for s in text.splitlines()) if line]
