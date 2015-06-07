
__accumulated_length__ = 0

def setup(context, meta):
    return meta

def cleanup(summary):
    summary['extra'] = { 'accumulated-length' : __accumulated_length__ }

    return summary

def process_feature(feature):
    global __accumulated_length__

    locus = feature['locus']

    length = abs(locus['end'] - locus['start'])

    if length < 10000:
        return None

    __accumulated_length__ += length

    feature['length'] = length

    if 'comment' in feature:
        del feature['comment']

    return feature
