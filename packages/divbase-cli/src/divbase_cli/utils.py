"""Collection of utility functions for divbase-cli package that haven't found a better home"""


def format_file_size(size_bytes: int | float | None) -> str:
    """
    Converts a file size in bytes to a human-readable format.

    Uses powers of 1000 so KB, MB, GB, TB and not 1024 KiB, MiB, GiB, TiB.
    """
    if size_bytes is None:
        return "N/A"
    if size_bytes == 0:
        return "0 B"
    power = 1000
    n = 0
    power_labels = {0: "", 1: "K", 2: "M", 3: "G", 4: "T"}
    while size_bytes >= power and n < len(power_labels) - 1:
        size_bytes /= power
        n += 1
    return f"{size_bytes:.2f} {power_labels[n]}B"
