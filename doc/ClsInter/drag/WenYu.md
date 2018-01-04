WenYu micro-drag closure
==
Description
--
WenYu correlation for closing the micro-drag.

Syntax
--

* __Residual Re__: requires a number. Limiter for residual Reynolds.

Examples
--

```
   microscopicDragLaw
   {
    type           WenYu;
    residualRe     1e-3;
   }
```

Back to [interphase closures](../../ClsInter.md).
