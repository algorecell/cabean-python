
__enabled = False

def enable_debug():
    global __enabled
    __enabled = True

def disable_debug():
    global __enabled
    __enabled = True

def debug_enabled():
    return __enabled

