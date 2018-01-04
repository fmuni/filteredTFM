Stress Closures
==

Description
--

Stress closures in __eulerianFilteredTFM__ are specified in [phasePropertiesDict](phasePropertiesDict.md)
for each phase. The dispersed phase requires three closures (meso-scale, frictional and micro-scale) while the continuous phase requires two closures (meso-scale, micro-scale).

Syntax
--
Each phase in [phasePropertiesDict](phasePropertiesDict.md) requires a sub-dictionary named _stressClosure_ which contains sub-dictionaries for each _sub-closure_. We refer to the functional form of a specific closure (i.e., meso-scale, frictional or microScale) as _sub-Closure_.

Each sub-closure will require the follwing entry:
* __type__: requires a string. Name of the model to use (entring __none__ will disable the sub-closure).

In order to stabilize some closures, it may be necessary to enforce a smooth transition to the
fully dilute regime (i.e., without dispersed phase) in some regions. By default, the smoother activates where
the particle volume fraction is smaller than 5%. This threshold can be adjusted by including the following:

* __dispersedPhaseThreshold__: requires a number. It sets the smoothing threshold.


Examples
--
_stressClosure_ for particle phase (dispersed).
```
...
stressClosure
{
      mesoScale
      {
          type  SarkarMeso;
          phase dispersed;

          dispersedPhaseThreshold 0.001;
      }

      frictional
      {
          type SchneiderbauerFrictional;

          I0      0.297;
          muSt    0.3819;
          muC     0.6435;
      }

      microScale
      {
          type SarkarMicro;
      }
}
...

```

Available stress sub-closures
--




| micro-scale               | frictional | meso-scale |
|:---                       |:---        |:---        |
| [molecular](ClsStress/micro/ms.md)  | [CloeteFrictional](ClsStress/frictional/Cloete.md) | [SarkarMeso](ClsStress/meso/Sarkar.md) |
| [SarkarMicro](ClsStress/micro/Sarkar.md)  | [SchneiderbauerFrictional](ClsStress/frictional/Schneiderbauer.md) |

Back to [main index](01_main.md).
