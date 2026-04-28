
# Bernstien.py

**Calculates Bernstein basis values**
<br>- Input P (order)
<br>- Input t

**Outputs Bernstein basis values and the sum**

# Bezier.py
## Draw Bezier curve based on the input given

**Input**
```
x = [x_1,y_1],
    [x_2,y_2],
    [x_3,y_3],
       .  . 
       .  . 
       .  . 
    [x_p,y_p]

```

# Multinomial_C1.py 

# Multinomial_C2.py

# Derivative_Bernstein.py


<br>**Inputs the order and the derivatives, $d$**
<br>**Outputs the start(t=0) and end(t=1) derivative values**

### Example for $ p = 2 $

$ B_0^2 = (1 - t)^2 $

$ B_1^2 = 2t(1 - t) $

$ B_2^2 = t^2 $

#### Differentiate:

$ \frac{d}{dt} B_0^2 = -2(1 - t) $

$ \frac{d}{dt} B_1^2 = 2 - 4t $

$ \frac{d}{dt} B_2^2 = 2t $

#### Using the formula:

$ \frac{d}{dt} B_1^2 = 2 \left( B_0^1 - B_1^1 \right) $


