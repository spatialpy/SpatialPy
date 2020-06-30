import logging

# Setup logging
log = logging.getLogger()

def init_log(name, modelname):
    """Initialize logger.
    Args:
        name (str) of caller class
        modelname (str) used as filename
    """

    log.name = name
    fh = logging.FileHandler(filename=modelname + '.log', encoding="utf-8")
    ch = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    log.setLevel(logging.WARNING)
    ch.setLevel(logging.WARNING)
    log.addHandler(fh)
    log.addHandler(ch)
