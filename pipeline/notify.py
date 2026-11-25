"""Stage 9 (tail): hand a one-line heads-up to whatever notifier the operator configured.

The reports are the substance; this just says "a scan finished, go look." Nothing is sent
unless `NOTIFY_CMD` is set. When it is, the command is run with the message appended as its
final argument and never through a shell, so the message cannot be interpreted:

    NOTIFY_CMD="ntfy publish genome"            # → ntfy publish genome "<message>"
    NOTIFY_CMD="/usr/local/bin/my-notifier"     # → my-notifier "<message>"

Failure to notify never fails the pipeline.
"""

from __future__ import annotations

import shlex
import subprocess

from . import config
from .util import log


def notify(message: str) -> bool:
    """Run `NOTIFY_CMD <message>`. True if it ran and exited 0; False (never raises) otherwise,
    including when no notifier is configured."""
    if not config.NOTIFY_CMD:
        return False
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
