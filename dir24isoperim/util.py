class Log: lvl = 0

def log(msg, indent=-1, end="\n"):
    if indent < 0: indent = Log.lvl
    print(indent*3*" " + msg, end=end)

def err(msg): log(FMT_FAIL%"Error: " + msg)

def warn(msg): log(FMT_WARN%"Warning: " + msg)

FMT_PASS = "\033[1;92m%s\033[0m"
FMT_FAIL = "\033[1;91m%s\033[0m"
FMT_WARN = "\033[1;93m%s\033[0m"

