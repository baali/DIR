Overview
--------

Idea is to search document images for particular image patch(of a
word) without doing OCR. Available literature for techniques which I
am following:

- `Keyword spotting in document images through word shape coding http://www.comp.nus.edu.sg/~tancl/publications/c2009/icdar09-baishuyong.pdf`_
- `Searching in document images cvit.iiit.ac.in/papers/jawahar04search.pdf`_

  
I tried to use popular feature extraction and matching algorithms like
SURF but didn't get good results yet. I have to explore more on how to
use them more efficiently and confirm if they can be used or not.

Dependencies
------------

Compiling instructions are part of comments in source files. There is
one file using OpenCV surf features but other programs depends on
`leptonica http://www.leptonica.com/`_ for all image related
operations.
  
