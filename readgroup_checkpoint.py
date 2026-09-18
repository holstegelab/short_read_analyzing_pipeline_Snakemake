"""Share the canonical readgroup checkpoint across Snakemake modules.

Snakemake deliberately gives each workflow module a separate ``checkpoints``
namespace. The readgroup-dependent Stats and Kraken input functions need the
checkpoint registered by Aligner, not their own empty namespaces. Ordinary
Python modules are imported once per controller process, making this a small,
explicit bridge without importing Aligner's rules more than once.
"""

_checkpoint = None


def register(checkpoint):
    """Register the checkpoint proxy created while importing Aligner.smk."""
    global _checkpoint
    _checkpoint = checkpoint


def output(sample):
    """Return the sample manifest and preserve Snakemake reevaluation."""
    if _checkpoint is None:
        raise RuntimeError(
            "get_readgroups checkpoint has not been registered; import the "
            "Aligner rules before readgroup-dependent modules"
        )
    return _checkpoint.get(sample=sample).output[0]
