#!/bin/bash
#
# Call from shell: tools/format.sh /path/to/file.cpp /path/to/some/other/file.cpp
# Call from vim for current buffer: :%!tools/format.sh
astyle \
    -s4     `# 4 space indent` \
    -S      `# Indent case: inside switch` \
    -O      `# Don't break one-line blocks` \
    -o      `# Don't break complex one-line expressions` \
    -k1     `# Attach pointer/reference modifier to type, not name (i.e. char** argv)` \
    -H      `# Insert space padding after if/for/while etc` \
    -p      `# Insert padding around operators` \
    -U      `# No padding inside parens` \
    -n      `# No backup file` \
    "$@"
