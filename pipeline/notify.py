"""Stage 9 (tail): hand a one-line heads-up to whatever notifier the operator configured.

The reports are the substance; this just says "a scan finished, go look." Nothing is sent
unless `NOTIFY_CMD` is set. When it is, the command is run with the message appended as its
final argument and never through a shell, so the message cannot be interpreted:

    NOTIFY_CMD="ntfy publish genome"            # → ntfy publish genome "<message>"
    NOTIFY_CMD="/usr/local/bin/my-notifier"     # → my-notifier "<message>"

Failure to notify never fails the pipeline.

**Quiet hours.** With `NOTIFY_QUIET_HOURS=00:00-07:00` a message raised in that window (local time
in `NOTIFY_TZ`, default Europe/London) is appended to `$GENOMES_ROOT/.notify_queue.jsonl` instead
of being sent. The queue goes out, oldest first, ahead of the next message sent outside the window,
or when `run.py notify-flush` runs (a timer just after the window ends). A message that fails to
send stays queued for the next attempt.
"""

from __future__ import annotations

import datetime as _dt
import json
import os
import shlex
import subprocess

from . import config
from .util import log


def _send(message: str) -> bool:
    try:
        argv = shlex.split(config.NOTIFY_CMD) + [message]
        r = subprocess.run(argv, capture_output=True, text=True, timeout=config.NOTIFY_TIMEOUT)
    except Exception as e:  # noqa: BLE001 — notification is best-effort
        log.warning("Notification failed (%s): %s", e, message)
        return False
    if r.returncode != 0:
        log.warning("Notification command exited %d (%s): %s",
                    r.returncode, (r.stderr or "").strip()[:200], message)
        return False
    log.info("Notified: %s", message)
    return True


def _now() -> _dt.datetime:
    try:
        from zoneinfo import ZoneInfo
        return _dt.datetime.now(ZoneInfo(config.NOTIFY_TZ))
    except Exception:  # noqa: BLE001 — no tz database: local time is the best guess
        return _dt.datetime.now()


def quiet_now(now: _dt.datetime | None = None) -> bool:
    """Is `now` (default: the current time in NOTIFY_TZ) inside NOTIFY_QUIET_HOURS?"""
    spec = (config.NOTIFY_QUIET_HOURS or "").strip()
    if not spec:
        return False
    try:
        a, b = (_dt.time.fromisoformat(x.strip()) for x in spec.split("-", 1))
    except ValueError:
        log.warning("NOTIFY_QUIET_HOURS=%r is not HH:MM-HH:MM; ignoring it.", spec)
        return False
    t = (now or _now()).time()
    return (a <= t < b) if a <= b else (t >= a or t < b)


def _queue_path():
    return config.GENOMES_ROOT / ".notify_queue.jsonl"


def _queued() -> list[dict]:
    p = _queue_path()
    out = []
    try:
        for line in p.read_text().splitlines():
            try:
                out.append(json.loads(line))
            except ValueError:
                continue
    except OSError:
        pass
    return out


def _write_queue(items: list[dict]) -> None:
    p = _queue_path()
    if not items:
        p.unlink(missing_ok=True)
        return
    tmp = p.with_name(f"{p.name}.{os.getpid()}.part")
    tmp.write_text("".join(json.dumps(i) + "\n" for i in items))
    os.replace(tmp, p)


def flush() -> int:
    """Send what quiet hours held back, oldest first; returns how many went out. Anything that
    fails to send stays queued. Does nothing inside quiet hours."""
    if not config.NOTIFY_CMD or quiet_now():
        return 0
    items, kept, sent = _queued(), [], 0
    for item in items:
        msg = item.get("message", "")
        when = item.get("queued", "")
        if msg and _send(f"{msg} (held overnight, raised {when})" if when else msg):
            sent += 1
        elif msg:
            kept.append(item)
    _write_queue(kept)
    return sent


def notify(message: str) -> bool:
    """Run `NOTIFY_CMD <message>`, or queue it inside quiet hours. True if it was sent (exit 0)
    or queued; False (never raises) otherwise, including when no notifier is configured."""
    if not config.NOTIFY_CMD:
        return False
    if quiet_now():
        items = _queued()
        items.append({"message": message, "queued": _now().strftime("%Y-%m-%d %H:%M")})
        try:
            config.GENOMES_ROOT.mkdir(parents=True, exist_ok=True)
            _write_queue(items)
            log.info("Quiet hours (%s %s): notification queued: %s",
                     config.NOTIFY_QUIET_HOURS, config.NOTIFY_TZ, message)
            return True
        except OSError as e:
            log.warning("Could not queue notification (%s): %s", e, message)
            return False
    flush()
    return _send(message)
