import json

_JSON_CHARACTER_REPLACEMENT_MAPPING = [

    ('<', '\\u%04x' % ord('<')),
    ('>', '\\u%04x' % ord('>')),
    ('&', '\\u%04x' % ord('&')),
]

class _JsonEncoderForHtml(json.JSONEncoder):
  def encode(self, o):
    chunks = self.iterencode(o, _one_shot=True)
    if not isinstance(chunks, (list, tuple)):
      chunks = list(chunks)
    return ''.join(chunks)

  def iterencode(self, o, _one_shot=False):
    chunks = super(_JsonEncoderForHtml, self).iterencode(o, _one_shot)
    for chunk in chunks:
      for (character, replacement) in _JSON_CHARACTER_REPLACEMENT_MAPPING:
        chunk = chunk.replace(character, replacement)
      yield chunk

def dump(*args, **kwargs):
    kwargs['cls'] = _JsonEncoderForHtml
    return json.dump(*args, **kwargs)

def dumps(*args, **kwargs):
    kwargs['cls'] = _JsonEncoderForHtml
    return json.dumps(*args, **kwargs)
