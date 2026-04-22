from typing import Literal, TypeAlias

QuoteChar: TypeAlias = Literal["'", '"']


def format_file_size(size_bytes: int | float | None, decimals: int = 2) -> str:
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
    return f"{size_bytes:.{decimals}f} {power_labels[n]}B"


def split_semicolon_bcftools_command_segments(command: str) -> list[str]:
    """
    Split a user provided bcftools command string on semicolons while respecting quoted substrings.

    bcftools view allows semicolons for certain options (e.g. -i 'FILTER="A;B"'), so we need to ensure that we only
    split on semicolons that are not inside quotes.

    --command "view -i 'FILTER=\"A;B\"'; view -r 1:1-1000"

    Slides character by character through the command string to identify the true semicolon delimiters, whilst ignoring semicolons inside quotes.
    The command strings are not that long, so we can afford to do this on a character level.
    """
    command_segments: list[str] = []
    segment_chars: list[str] = []

    single_quote = "'"
    double_quote = '"'

    # Possible states None = outside quotes, "'" = inside single-quoted substring, '"' = inside double-quoted substring.
    open_quote_char: QuoteChar | None = None
    is_escaped = False

    for char in command:
        if is_escaped:
            segment_chars.append(char)
            is_escaped = False
            continue

        # Handle backslash-escaped characters in the received command text.
        # In single-quoted regions, backslashes are treated as literal chars.
        if char == "\\" and open_quote_char != single_quote:
            segment_chars.append(char)
            is_escaped = True
            continue

        if char in (double_quote, single_quote):
            segment_chars.append(char)
            if open_quote_char is None:
                open_quote_char = char
            elif open_quote_char == char:
                open_quote_char = None
            continue

        if char == ";" and open_quote_char is None:
            command_segments.append("".join(segment_chars))
            segment_chars = []
            continue

        segment_chars.append(char)

    command_segments.append("".join(segment_chars))
    return command_segments
