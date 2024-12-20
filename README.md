# Constructing a Boolean model of HER2+ breast cancer

Human epidermal growth factor receptor 2 (HER2)-positive breast cancer is characterised by overexpression of the HER2 receptor or its downstream kinases, resulting in an aggressive phenotype that impacts 25% of breast cancer patients. HER2 is localised to the plasma membrane in several tissues and belongs to multiple cell-surface receptor complexes, making it an easily exploitable molecular target. Consequently, a preferred treatment for HER2+ breast cancer is the humanized IgG1 monoclonal antibody, trastuzumab, that targets HER2. Trastuzumab has contributed to patient survival rates increasing to over 90%, however, HER2+ breast cancer shortly became susceptible to trastuzumab resistance with frequent tumour relapses. The partial cause of the trastuzumab resistance is the activity of several dual-specificity phosphatases (DUSPs). 

In this repository, I'll try to recreate a Boolean model constructed by [Buiga and colleagues](https://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-018-0534-5) (2018), and simulate *in vitro* HER2+ breast cancer cell survival outcomes when treated with trastuzumab, and DUSP16 inhibitor. Several DUSPs are already established as inducers of trastuzumab resistance, including DUSP4, DUSP6 and DUSP8, but not as much is known about DUSP16, despite some evidence showing similar cellular responses in the context of HER2+ breast cancer, which is why I wanted to construct my model focusing on DUSP16.

## 🧬 Regulatory mechanisms

Based on the Buiga et a., (2018) paper and some of their subsequent work, the following regulatory mechanisms were identified:
- Trastuzumab inhibits HER2
- HER2 activates ERK1/2
- HER2 activates JNK1/2
- HER2 activates P38
- ERK1/2 activates DUSP16
- DUSP16 inhibits JNK1/2
- DUSP16 inhibits P38
- ERK1/2 activates Survival
- JNK1/2 and P38 together inhibit Survival

This can be rewritten using regulatory logic rules, with each protein, pharmacological compound or the outcome (e.g. Survival) existing as individual nodes, as below. 
Node | Regulatory logic rule
--- | ---
HER2 | HER2<sup>t+1</sup> = not Trastuzumab<sup>t</sup>
ERK1/2 | ERK1/2<sup>t+1</sup> = HER2<sup>t</sup>
JNK1/2 | JNK1/2<sup>t+1</sup> = HER2<sup>t</sup> or not DUSP16<sup>t</sup>
P38 | P38<sup>t+1</sup> = HER2<sup>t</sup> or not DUSP16<sup>t</sup>
DUSP | DUSP16<sup>t+1</sup> = ERK1/2<sup>t</sup> or not DUSP16_inhibitor<sup>t</sup>
Survival | Survival<sup>t+1</sup> = ERK1/2<sup>t</sup> or not (JNK1/2<sup>t</sup> and P38<sup>t</sup>)

To help visualise this Boolean network, I used yEd Live to produce the following plot, where arrow strokes indicate activation, tick mark strokes indicate inhibition and an orange star 
indicates cumulative action of P38 and JNK1/2.<br/>

<img src="https://github.com/user-attachments/assets/370e36a9-1531-4f0c-b8e5-7c98e27ee8cc" width="425" />

## Boolean model development

To develop this model, we are going to load the following packages:
```
import boolean2 as b2
import boolean2, pylab
from boolean2 import util
import matplotlib.pyplot as plt
```

Note: boolean2 package (belonging to BooleanNet) can only be run using Python 2, so you will have to ensure you are using Python 2 kernel instead of Python 3 kernel. 

Considering the model rules outlined above, we are going to define our HER2+ breast cancer model (prior to any treatment) as follows:
```
print("Prior to Trastuzumab treatment when HER2 is ON:")

model_definition = """

# Initial values 
HER2 = True
ERK12 = Random
JNK12 = Random
P38 = Random
DUSP16 = Random
Survival = True

# Update rules
ERK12* = HER2
JNK12* = HER2 or not DUSP16
P38* = HER2 or not DUSP16
DUSP16* = ERK12
Survival* = ERK12 or not (JNK12 and P38)
"""

model = b2.Model(text = model_definition, mode = 'async')
model.initialize()
model.iterate(steps = 100) # Change to 10 when exporting the figures
```

This iterated a model simulation 100 times, and we can now assess whether the model converges to a steady state or a limit cycle using the following code:
```
for state in model.states:
    print state.HER2, state.ERK12, state.JNK12, state.P38, state.DUSP16, state.Survival
    
model.report_cycles()
```

This will print the following output, informing us of the state of each node (i.e. is it active (True) or inactive (False)), and that the model converges to a steady state at index 1.

>Prior to Trastuzumab treatment when HER2 is ON:<br/>
>True True True True False True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>...<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>True True True True True True<br/>
>Steady state starting at index 1 -> State: DUSP16=True, ERK12=True, HER2=True, JNK12=True, P38=True, Survival=True

This exact output can be visualised, and to do that, we will first change the number of iterations to 10 to ensure we get a tidy plot, and use the following code to plot a simulation heatmap:
```
image = list()
for node in model.data:
    image.append(model.data[node])
plt.yticks(range(0,6), model.data)
plt.imshow(image, cmap=plt.cm.get_cmap('RdYlGn'), interpolation='none')

plt.savefig('Fig1-no-inhibitors.png')
```

Now we have generated the following heatmap. The green regions represent the active or true ("ON") states of the nodes on the Y-axis, whereas the red regions represent inactive or false ("OFF") states, over a course of 10 iterations, considered arbitrary time points.<br/>

![Fig1-no-inhibitors](https://github.com/user-attachments/assets/4d1d28e3-37d6-4668-9ab6-11590cf72740)

From this model simulation, we can see that the HER2+ breast cancer cell survival rate is quite high, as the Survival node is active throughout all time points, representing a 100% survival rate. 

## Trastuzumab treatment simulation on the developed HER2+ breast cancer model

We will now alter the model by adding trastuzumab to the simulation and changing the initial values of certain nodes. Specifically, we are going to set the initial value of HER2 node to False because of the addition of trastuzumab and the initial value of the Survival node to random. 
```
print("After Trastuzumab treatment when HER2 is OFF:")

model_definition = """

# Initial values 
Trastuzumab = True
HER2 = False
ERK12 = Random
JNK12 = Random
P38 = Random
DUSP16 = Random
Survival = Random

# Update rules
HER2* = not Trastuzumab
ERK12* = HER2
JNK12* = HER2 or not DUSP16
P38* = HER2 or not DUSP16
DUSP16* = ERK12
Survival* = ERK12 or not (JNK12 and P38)
"""

model = b2.Model(text = model_definition, mode = 'async')
model.initialize()
model.iterate(steps = 100) # Change to 10 when exporting figures

# Assess whether the model converges to a steady state or a limit cycle
for state in model.states:
    print state.Trastuzumab, state.HER2, state.ERK12, state.JNK12, state.P38, state.DUSP16, state.Survival
    
model.report_cycles()
```

The output of this simulation informs us that the altered model converges to a steady state at index 3, and can be visualised using the same code as above:
```
image = list()
for node in model.data:
    image.append(model.data[node])
plt.yticks(range(0,7), model.data)
plt.imshow(image, cmap=plt.cm.get_cmap('RdYlGn'), interpolation='none')

plt.savefig('Fig2-HER2-inhibition.png')
```

Using this code chunk, we generate the following heatmap.

![Fig2-HER2-inhibition](https://github.com/user-attachments/assets/86664a38-08d8-48cf-88cc-570ad79205ef)

Looking at this plot, we notice a significant difference compared to the initial plot that models cell survival without any treatment. By treating the cells with trastuzumab, the model predicts that cells' survival will only occur within the initial three time points, followed by apoptosis, decreasing the overall HER2+ breast cancer cell survival rate by 72.72%. 

## Simulation of the combined Trastuzumab + DUSP16 inhibitor therapy
To assess how the survival outcomes change with the addition of a DUSP16 inhibitor, we will alter our initial values and model rules once again in the code below. 

```
print("After DUSP inhibition:")

model_definition = """

# Initial values 
Trastuzumab = True
DUSP_inhibitor = True
HER2 = False
ERK12 = Random
JNK12 = Random
P38 = Random
DUSP16 = False
Survival = Random

# Update rules
HER2* = not Trastuzumab
ERK12* = HER2
DUSP16* = not DUSP_inhibitor
JNK12* = HER2 or not DUSP16
P38* = HER2 or not DUSP16
Survival* = ERK12 or not (JNK12 and P38)
"""

model = b2.Model(text = model_definition, mode = 'async')
model.initialize()
model.iterate(steps = 10) # Change to 10 when exporting figures

# Assess whether the model converges to a steady state or a limit cycle
for state in model.states:
    print state.Trastuzumab, state.DUSP_inhibitor, state.HER2, state.ERK12, state.JNK12, state.P38, state.DUSP16, state.Survival
    
model.report_cycles()
```

This simulation results in the following (visual) output:

![Fig3-double-inhibition](https://github.com/user-attachments/assets/92326717-4f69-4426-96d7-63ee6028817f)

From this output, we notice that the combination of trastuzumab and DUSP16 inhibitor converge the model towards a steady state starting at index 1, and results in the Survival node being active at 1 out of 11 time points, bringing the HER2+ breast cancer cell survival rate further down to 9.09%. Thus, the model predicts that combining trastuzumab and DUSP16 inhibitor will result in a 90.90% apoptosis rate.

Boolean models like the ones developed in this repository, make a useful tool to predict qualitative states and interactions of molecular components within the biological process of interest when quantitative data describing the said biological process is lacking. Using stochastic models with asynchronous updates, as we did here, is frequently reported in research articles as a way of simulating apoptotic signalling and drug resistance in cancer cells. However, it is important to mention that network heatmaps generated using this type of modelling do not represent a fixed state due to the asynchronous nature of model updates. This means that every time a model is simulated, the output will have a slightly different prediction. Although rarely significant, this difference in output predictions implies that they are not definitive but only approximate. 



