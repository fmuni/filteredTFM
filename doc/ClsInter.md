Interphase closures
==

Description
--

Interphase closures are used to represent interaction between different phase fields.
Specifically momentum and heat transfer (drag and heat transfer coefficient).

Syntax
--

They are specified in [phasePropertiesDict](phasePropertiesDict.md) as:

* __drag__: sub-dictionay for the drag coefficient.
 * __microscopicDragLaw__: sub-dictionary for the microscopic closure.
 * __heterogeneousCorrection__: sub-dictionary for the heterogeneous correction.


* __heatTransfer__: sub-dictionay for the heat transfer coefficient.
 * __microscopicNusselt__: sub-dictionary for the microscopic closure.
 * __heterogeneousCorrection__: sub-dictionary for the heterogeneous correction.

Examples
--
```
drag
{
       microscopicDragLaw
       {
        type            GidaspowErgunWenYu;
        residualRe      1e-3;
       }

       heterogeneousCorrection
       {
        type        SchneiderbauerPirker;
       }

}


heatTransfer
{
      microscopicNusselt
      {
        type            RanzMarshall;
        residualAlpha   1e-3;
      }

      heterogeneousCorrection
      {
       type        none;
      }

}
```
Available interphase closures

| drag microscopic | drag heterogeneous correction | Nusselt microscopic | Nusselt heterogeneous correction |
|:-- |:-- |:-- |:-- |
| [Ergun](ClsInter/drag/Ergun.md)  | [homogeneous](ClsInter/hdrag/homo.md)  | [RanzMarshall](ClsInter/heat/RM.md)| [homogeneous](ClsInter/hdrag/homo.md) |
| [Gibilaro](ClsInter/drag/Gibilaro.md)  | [Cloete](ClsInter/hdrag/Cloete.md)  | [Spherical](ClsInter/heat/Spherical.md) |
| [GidaspowErgunWenYu](ClsInter/drag/GidaspowErgunWenYu.md)  | [Igci](ClsInter/hdrag/Igci.md) ||
| [GidaspowSchillerNaumann](ClsInter/drag/GidaspowSchillerNaumann.md)  | [Milioli](ClsInter/hdrag/Milioli.md) ||
| [IshiiZuber](ClsInter/drag/IshiiZuber.md)  | [SchneiderbauerPirker](ClsInter/hdrag/SP.md)  ||
| [Lain](ClsInter/drag/Lain.md)  |  ||
| [SchillerNaumann](ClsInter/drag/SchillerNaumann.md)  |||
| [SyamlalOBrien](ClsInter/drag/SyamlalOBrien.md)  |||
| [WenYu](ClsInter/drag/WenYu)  |||

Back to [main index](01_main.md).
