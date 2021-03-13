#!/usr/bin/env python3

import sys
import re

def formatKey(key):
    if key.startswith("Example:"):
        return key
    else:
        return key + ":"

key = ""
isMarkdown = False

for line in sys.stdin:
    line = line.rstrip()

    tag = line.strip()
    if re.search("^<(\w|/)", tag) != None:
        if tag.startswith("<div class=\"panel-heading\">"):
            print("### " + re.match("^<.*?>(.*)</div>$", tag).group(1))
        elif tag.startswith("</") or tag.startswith("<div ") or tag.startswith("<dl"):
            if tag == "</dd>" and isMarkdown:
                print("</pre>")
                isMarkdown = False
            pass
        elif tag.startswith("<dt>"):
            key = re.match("^<.*>(.*?)<.*>$", tag).group(1)
        elif tag == "<dd markdown=\"1\">":
            print("* **" + formatKey(key) + "**")
            print("<pre>")
            isMarkdown = True
        elif tag.startswith("<dd>"):
            print("* **" + formatKey(key) + "** " + re.match("^<.*?>(.*?)(?:</dd>)?$", tag).group(1))
        elif tag.startswith("<dd "):
            print("* **" + formatKey(key) + "**")
        else:
            print(line)
    else:
        print(line)

