# Constructing a Boolean model of HER2+ breast cancer

Human epidermal growth factor receptor 2 (HER2)-positive breast cancer is characterised by overexpression of the HER2 receptor or its downstream kinases, resulting in an aggressive phenotype that impacts 25% of breast cancer patients. HER2 is localised to the plasma membrane in several tissues and belongs to multiple cell-surface receptor complexes, making it an easily exploitable molecular target. Consequently, a preferred treatment for HER2+ breast cancer is the humanized IgG1 monoclonal antibody, trastuzumab, that targets HER2. Trastuzumab has contributed to patient survival rates increasing to over 90%, however, HER2+ breast cancer shortly became susceptible to trastuzumab resistance with frequent tumour relapses. The partial cause of the trastuzumab resistance is the activity of several dual-specificity phosphatases (DUSPs). 

In this repository, I'll try to recreate a Boolean model constructed by [Buiga and colleagues](https://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-018-0534-5) (2018), and simulate *in vitro* HER2+ breast cancer cell survival outcomes when treated with trastuzumab, and DUSP16 inhibitor. Several DUSPs are already established as inducers of trastuzumab resistance, including DUSP4, DUSP6 and DUSP8, but not as much is known about DUSP16, despite some evidence showing similar cellular responses in the context of HER2+ breast cancer, which is why I wanted to construct my model focusing on DUSP16.

## 🧬 Regulatory mechanisms

Based on the Buiga et a., (2018) paper and some of their subsequent work, the following regulatory mechanisms were identified:
- Herceptin inhibits HER2
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
HER2 | HER2<sup>t+1</sup> = not Herceptin<sup>t</sup>
ERK1/2 | ERK1/2<sup>t+1</sup> = HER2<sup>t</sup>
JNK1/2 | JNK1/2<sup>t+1</sup> = HER2<sup>t</sup> or not DUSP16<sup>t</sup>
P38 | P38<sup>t+1</sup> = HER2<sup>t</sup> or not DUSP16<sup>t</sup>
DUSP | DUSP16<sup>t+1</sup> = ERK1/2<sup>t</sup> or not DUSP16_inhibitor<sup>t</sup>
Survival | Survival<sup>t+1</sup> = ERK1/2<sup>t</sup> or not (JNK1/2<sup>t</sup> and P38<sup>t</sup>)

To help visualise this Boolean network, I used yEd Live to produce the following plot.<br/>

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
print("Prior to Herceptin treatment when HER2 is ON:")

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

To help visualise the simulations, we will first change the number of iterations to 10 to ensure we get a tidy plot, and use the following code to plot a simulation heatmap:
```
image = list()
for node in model.data:
    image.append(model.data[node])
plt.yticks(range(0,6), model.data)
plt.imshow(image, cmap=plt.cm.get_cmap('RdYlGn'), interpolation='none')

plt.savefig('Fig1-no-inhibitors.png')
```

This generates the following heatmap. The green regions represent the active ("ON") state of the nodes on the Y-axis, whereas the red regions represent inactive ("OFF") states, over a course of 10 iterations, considered as arbitrary time points.<br/>
![Fig1-no-inhibitors](https://github.com/user-attachments/assets/4d1d28e3-37d6-4668-9ab6-11590cf72740)




