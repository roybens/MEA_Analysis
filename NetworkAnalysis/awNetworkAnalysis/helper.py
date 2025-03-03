import sys
import functools

# Global variable to track indentation
_INDENT_MODE = False
_INDENT_LEVEL = 0

def _indented_print(*args, **kwargs):
    """Custom print function that respects indentation settings."""
    indent = " " * (_INDENT_LEVEL * 4)  # Each level = 4 spaces
    built_in_print(indent + " ".join(map(str, args)), **kwargs)

def indent_mode_on(level=1):
    """Enable indentation mode with a specified level."""
    global _INDENT_MODE, _INDENT_LEVEL, built_in_print
    if not _INDENT_MODE:
        built_in_print = print  # Store original print function
        sys.modules["builtins"].print = _indented_print  # Override print function
        _INDENT_MODE = True
    _INDENT_LEVEL = level

def indent_mode_off():
    """Disable indentation mode and restore original print behavior."""
    global _INDENT_MODE
    if _INDENT_MODE:
        sys.modules["builtins"].print = built_in_print  # Restore original print
        _INDENT_MODE = False
        
def indent_increase(increment=1):
    """Increase the current indentation level."""
    global _INDENT_LEVEL, _INDENT_MODE
    if _INDENT_MODE:
        _INDENT_LEVEL += increment
    else:
        # enable indentation mode with the specified increment
        indent_mode_on(level=increment)
        
def indent_decrease(decrement=1):
    """Decrease the current indentation level."""
    global _INDENT_LEVEL, _INDENT_MODE
    if _INDENT_MODE:
        _INDENT_LEVEL -= decrement
        if _INDENT_LEVEL < 0:
            _INDENT_LEVEL = 0
    else:
        # print warning
        print('Warning: Indentation mode is not enabled. Enable it using indent_mode_on() or indent_increase() functions.')