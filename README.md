### What is it?
This is a reimplementation of [Zopfli](https://github.com/google/zopfli) in C#. It's a straightforward port from C, without relying on many C# specific features.

### Why rewrite in C#?
Mainly because I wanted to teach myself C#.
Other reason was because of the old blog posts from [Raymon Chen](https://devblogs.microsoft.com/oldnewthing/20050510-55/?p=35673) and [Rico Mariani](https://docs.microsoft.com/en-us/archive/blogs/ricom/performance-quiz-6-chineseenglish-dictionary-reader) where the first implemented a dictionary in C++ and the latter reimplemented it in C#. What may be suprising, the C# version turned out to be not much slower than highly optimized C++ impementation. I wanted to perform similar experiment and zopfli is good candidate, because it's quite CPU heavy.

### How does it compare to original?
My version is 1.6 slower than zopfli compiled with Microsoft C/C++ compiler. Not bad in my opinion.

### Can it be quicker?
Maybe. I think it can be trivially parallelised. Additionally there is at least one optimization left, but my current knowledge of C# isn't good enough to do it. Maybe I manage it in the future.

### Should I use it in production?
No. Use original compiled with good compiler. It's still noticeably quicker.

## Chapter 2: The generative AI era

### What is it this time?
Since recently there is a lot of talk about vibe coding and agentic engineering I decided to have some fun with it. I dug out this project and let Claude Opus work on it. Goal was mainly optimization, but we got some additional perks as we went by.

### How does it compare to original now?
The new version is around 50% quicker. But this time I also used the newest clang version on the C original with these flags ``-O3 -march=native -mtune=native -flto -fuse-ld=lld -fomit-frame-pointer -DNDEBUG``. These should have produced quite optimized binary. Results? C# version is still around 1.6 slower than C.

### Should I use it in production?
Probably still not. It is slower, but on the other hand thanks to Claude it has some features absent from original:

* Compressed file is streamed, so on big files memory usage is much smaller.
* My (ok, Claude's) version writes to temp file first, which is deleted if program is interrupted.

So yeah, you probably should use something else anyway.

## Bonus
If for some reason you want to see how I worked with Claude and how Claude itself worked, here is a [transcript](https://gist.github.com/Tisit/b1f96b9734db6811a25c190ac013fae7) documenting most of the work done on this codebase.