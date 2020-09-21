## What is it?
This is a reimplementation of [Zopfli](https://github.com/google/zopfli) in C#. It's a straightforward port from C, without relying on many C# specific features.

## Why rewrite in C#?
Mainly because I wanted to teach myself C#.
Other reason was because of the old blog posts from [Raymon Chen](https://devblogs.microsoft.com/oldnewthing/20050510-55/?p=35673) and [Rico Mariani](https://docs.microsoft.com/en-us/archive/blogs/ricom/performance-quiz-6-chineseenglish-dictionary-reader) where the first implemented a dictionary in C++ and the latter reimplemented it in C#. What may be suprising, the C# version turned out to be not much slower than highly optimized C++ impementation. I wanted to perform similar experiment and zopfli is good candidate, because it's quite CPU heavy.

## How does it compare to original?
My version is 1.6 slower than zopfli compiled with Microsoft C/C++ compiler. Not bad in my opinion.

## Can it be quicker?
Maybe. I think it can be trivially parallelised. Additionally there is at least one optimization left, but my current knowledge of C# isn't good enough to do it. Maybe I manage it in the future.

## Should I use it in production?
No. Use original compiled with good compiler. It's still noticeably quicker.
