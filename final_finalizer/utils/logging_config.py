#!/usr/bin/env python3
"""
Logging configuration for final_finalizer.

Provides a consistent logging setup with custom formatting that matches
the existing [info], [warn], [error], [done] pattern with ANSI formatting
for phase announcements.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional


# ANSI escape codes for formatting
BOLD = "\033[1m"
RESET = "\033[0m"
# Colors for different log levels (used in terminal output)
GREEN = "\033[32m"
YELLOW = "\033[33m"
RED = "\033[31m"
CYAN = "\033[36m"


class FinalFinalizerFormatter(logging.Formatter):
    """Custom formatter matching [info], [warn], [error], [done] pattern."""

    def __init__(self, use_color: bool = True):
        super().__init__()
        self.use_color = use_color

    def format(self, record: logging.LogRecord) -> str:
        # Handle custom log levels
        level_name = record.levelname.lower()

        # Map log levels to tags
        tag_map = {
            "debug": "debug",
            "info": "info",
            "warning": "warn",
            "error": "error",
            "critical": "error",
            "phase": "info",  # Custom level for phase announcements
            "done": "done",   # Custom level for completion messages
        }

        tag = tag_map.get(level_name, level_name)

        # Format message
        msg = record.getMessage()

        # Check for phase or done messages (custom attribute)
        if hasattr(record, 'is_phase') and record.is_phase:
            if self.use_color:
                return f"{BOLD}{CYAN}[{tag}] {msg}{RESET}"
            return f"[{tag}] {msg}"
        elif hasattr(record, 'is_done') and record.is_done:
            if self.use_color:
                return f"{GREEN}[done]{RESET} {msg}"
            return f"[done] {msg}"
        else:
            return f"[{tag}] {msg}"


class StderrHandler(logging.StreamHandler):
    """Handler that always writes to stderr."""

    def __init__(self):
        super().__init__(sys.stderr)


# Custom log level for phase announcements
PHASE_LEVEL = 25  # Between INFO (20) and WARNING (30)
DONE_LEVEL = 21   # Just above INFO


def _add_custom_levels():
    """Add custom log levels for phase and done messages."""
    logging.addLevelName(PHASE_LEVEL, "PHASE")
    logging.addLevelName(DONE_LEVEL, "DONE")


def _phase(self, message, *args, **kwargs):
    """Log a phase announcement (bold formatting in terminal)."""
    if self.isEnabledFor(PHASE_LEVEL):
        # Add custom attribute for formatting
        kwargs.setdefault('extra', {})['is_phase'] = True
        self._log(PHASE_LEVEL, message, args, **kwargs)


def _done(self, message, *args, **kwargs):
    """Log a completion message."""
    if self.isEnabledFor(DONE_LEVEL):
        kwargs.setdefault('extra', {})['is_done'] = True
        self._log(DONE_LEVEL, message, args, **kwargs)


# Add custom methods to Logger class
logging.Logger.phase = _phase
logging.Logger.done = _done


def setup_logging(
    verbose: bool = False,
    quiet: bool = False,
    log_file: Optional[str] = None,
) -> None:
    """Configure logging for final_finalizer.

    Sets up logging with the custom formatter that matches the existing
    [info], [warn], [error], [done] output pattern.

    Args:
        verbose: Enable DEBUG level logging
        quiet: Suppress INFO messages (only show warnings and errors)
        log_file: Optional path to also write logs to a file
    """
    _add_custom_levels()

    # Determine log level
    if verbose:
        level = logging.DEBUG
    elif quiet:
        level = logging.WARNING
    else:
        level = logging.INFO

    # Get root logger for final_finalizer
    root_logger = logging.getLogger("final_finalizer")
    root_logger.setLevel(logging.DEBUG)  # Allow all messages through
    root_logger.handlers.clear()  # Remove any existing handlers

    # Console handler (stderr)
    use_color = sys.stderr.isatty()
    console_handler = StderrHandler()
    console_handler.setLevel(level)
    console_handler.setFormatter(FinalFinalizerFormatter(use_color=use_color))
    root_logger.addHandler(console_handler)

    # File handler (optional)
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_path)
        file_handler.setLevel(logging.DEBUG)  # Always log everything to file
        file_handler.setFormatter(FinalFinalizerFormatter(use_color=False))
        root_logger.addHandler(file_handler)


def get_logger(name: str) -> logging.Logger:
    """Get a logger for a specific module.

    Args:
        name: Module name (e.g., "cli", "classifier")

    Returns:
        Configured logger instance with custom phase() and done() methods
    """
    return logging.getLogger(f"final_finalizer.{name}")
