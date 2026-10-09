"""Restore the companion Supplementary Code HTML; no analysis code is run."""
from pathlib import Path, PurePosixPath
from html.parser import HTMLParser
import argparse, base64, gzip, hashlib

class PackageParser(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.files = []
        self.current = None
    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if tag == "pre" and "data-file" in attributes:
            if self.current is not None:
                raise ValueError("Nested source block")
            self.current = (attributes, [])
    def handle_data(self, value):
        if self.current is not None:
            self.current[1].append(value)
    def handle_endtag(self, tag):
        if tag == "pre" and self.current is not None:
            attrs, chunks = self.current
            text = "".join(chunks)
            if attrs["data-encoding"] in {"base64", "gzip-base64"}:
                data = base64.b64decode(text, validate=True)
                if attrs["data-encoding"] == "gzip-base64":
                    data = gzip.decompress(data)
            else:
                if attrs["data-newline"] == "crlf":
                    text = text.replace("\n", "\r\n")
                data = text.encode("utf-8")
                if attrs["data-bom"] == "1":
                    data = b"\xef\xbb\xbf" + data
            if hashlib.sha256(data).hexdigest() != attrs["data-sha256"]:
                raise ValueError("Checksum mismatch: " + attrs["data-file"])
            self.files.append((attrs["data-file"], data))
            self.current = None

def restore(source, destination):
    parser = PackageParser()
    parser.feed(Path(source).read_bytes().decode("utf-8"))
    parser.close()
    if not parser.files or parser.current is not None:
        raise ValueError("No complete package found")
    root = Path(destination).resolve()
    paths = []
    for name, data in parser.files:
        relative = PurePosixPath(name)
        if relative.is_absolute() or ".." in relative.parts or "\\" in name or ":" in name:
            raise ValueError("Unsafe path: " + name)
        target = root.joinpath(*relative.parts).resolve()
        if not target.is_relative_to(root) or target in paths:
            raise ValueError("Unsafe or duplicate path: " + name)
        if target.exists() and target.read_bytes() != data:
            raise FileExistsError("Refusing to overwrite a different file: " + name)
        paths.append(target)
    for target, (_, data) in zip(paths, parser.files):
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
    return len(paths)

if __name__ == "__main__":
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument("html")
    args.add_argument("output", help="New directory for the restored files")
    values = args.parse_args()
    count = restore(values.html, values.output)
    print(f"Restored {count} files; all SHA256 checks passed.")
