
# This variable will be used to accumulate the
# length (in basepairs) of all features that are
# not filtered out in `process_feature` below.
#
# Use this or a similar approach to share data
# across features.
__accumulated_length__ = 0

def setup(context, meta):
    # The context is always being output, but meta
    # could be omitted. How? It is quite simple:
    # returning (a possibly modified) meta Dict
    # will be output, but returning None will
    # output no meta information.
    return meta

def cleanup(summary):
    # Add the accumulated length of all features
    # to the summary:
    summary['my-user-data'] = { 'accumulated-length' : __accumulated_length__ }

    # The Dict `summary` needs to be returned in order
    # to appear in the output. If no summary should appear
    # in the output, then None needs to be returned.
    return summary

def process_feature(feature):
    global __accumulated_length__

    locus = feature['locus']

    length = abs(locus['end'] - locus['start']) + 1

    # If we want to filter out some features (remove them
    # from the output completely), then we can do so my
    # returning None.
    if length < 10000:
        return None

    __accumulated_length__ += length

    feature['my-calculated-length'] = length

    if 'comment' in feature:
        del feature['comment']

    return feature
