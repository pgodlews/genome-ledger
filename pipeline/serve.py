"""Local web server for the report dashboards.

Serves ONLY the report directories — raw genomes are never exposed.
stdlib only; bound to localhost by default — it is unauthenticated, so keep it LAN-only.
"""

from __future__ import annotations

from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse

from . import config
from .util import log

# The only web-servable directory is the HTML reports tree. Reports used to sit inside the
# stage directories, beside genotype data, which is why an allowlist of directories had to be
# paired with an allowlist of extensions; now nothing under reports/html/ is anything but a
# rendered page. raw/, aligned/, called/, duckdb/, the stage dirs and even reports/md/ are
# never reachable over HTTP.
_ALLOWED = ("reports/html",)
_TYPES = {".html": "text/html; charset=utf-8", ".md": "text/markdown; charset=utf-8",
          ".json": "application/json", ".txt": "text/plain; charset=utf-8",
          ".css": "text/css", ".woff2": "font/woff2"}


def content_type(target: Path) -> str | None:
    """Content-Type for a servable report file, or None if it must not be served.

    An extension allowlist, deliberately on top of the directory allowlist above. The
    servable tree holds only rendered pages today, but serving any file under an allowed dir
    (the old octet-stream fallback) is one misplaced file away from putting genotype data on
    an unauthenticated LAN port."""
    return _TYPES.get(target.suffix)


class _Handler(BaseHTTPRequestHandler):
    def log_message(self, *a):  # noqa: D401 — silence default stderr access log
        pass

    def _send(self, code: int, body: bytes, ctype: str = "application/json") -> None:
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        parsed = urlparse(self.path)
        rel = parsed.path.lstrip("/")
        if rel in ("", "index.html"):
            # One landing page, the static one `render` writes — the same file whether it
            # is opened from disk or reached here. Redirect rather than serve its bytes at
            # "/", because its links are relative to its own directory.
            if (config.REPORTS_DIR / "html" / "index.html").exists():
                self.send_response(302)
                self.send_header("Location", "/reports/html/index.html")
                self.send_header("Content-Length", "0")
                self.end_headers()
                return None
            return self._send(200, _NO_REPORTS.encode(), "text/html; charset=utf-8")
        root = config.GENOMES_ROOT.resolve()
        # Bookmarks from before the pages moved under reports/html/: /reports/<snapshot>/…
        # is answered with a redirect when that page exists at its new address.
        if rel.startswith("reports/") and not rel.startswith("reports/html/"):
            moved = (root / "reports" / "html" / rel[len("reports/"):]).resolve()
            html_root = (root / "reports" / "html").resolve()
            if html_root in moved.parents and moved.is_file():
                self.send_response(302)
                self.send_header("Location", "/reports/html/" + rel[len("reports/"):])
                self.send_header("Content-Length", "0")
                self.end_headers()
                return None
        target = (root / rel).resolve()
        allowed = [(root / a).resolve() for a in _ALLOWED]
        if not any(ar == target or ar in target.parents for ar in allowed):
            return self._send(403, b'{"error":"forbidden"}')   # blocks .. and raw dirs
        if not target.is_file():
            return self._send(404, b'{"error":"not found"}')
        ctype = content_type(target)
        if ctype is None:
            return self._send(403, b'{"error":"forbidden"}')
        return self._send(200, target.read_bytes(), ctype)

_NO_REPORTS = ("<!DOCTYPE html><meta charset=utf-8><title>Genome reports</title>"
               "<h1>Genome reports</h1><p>No reports yet. Run <code>scan</code> (or "
               "<code>render</code>) to build them.</p>")


def serve(host: str = "127.0.0.1", port: int = 8765) -> None:
    httpd = ThreadingHTTPServer((host, port), _Handler)
    log.info("Serving reports at http://%s:%d/  (Ctrl-C to stop)", host, port)
    if host not in ("127.0.0.1", "localhost"):
        log.warning("Bound to %s — this server is UNAUTHENTICATED; keep it LAN-only.", host)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        log.info("Stopped.")
        httpd.shutdown()
