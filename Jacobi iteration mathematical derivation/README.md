# Direct Computation of the nth Iteration in the Two-Variable Jacobi Method

**Field:** Numerical Analysis  
**Topic:** Iterative Methods and Recurrence Relations  
**Approach:** Direct computation using scalar recurrence relations

---

## 1. Introduction

The Jacobi method is an iterative numerical technique used to solve systems of simultaneous linear equations.

In its conventional implementation, successive approximations are calculated one iteration at a time. For example, calculating the fifth iteration normally requires calculating the first, second, third, and fourth iterations.

This project explores an algebraic shortcut for two-variable Jacobi systems.

By eliminating one variable and deriving a scalar recurrence relation, selected iterations can be calculated directly without computing every intermediate iteration.

The approach combines:

- Algebraic substitution.
- Scalar recurrence relations.
- Geometric progression.
- Closed-form mathematical expressions.

## 2. Mathematical Foundation

Consider the system:

$$
a_1x+b_1y=c_1
$$

$$
a_2x+b_2y=c_2
$$

Rearranging into Jacobi iteration form:

$$
x^{(k+1)}=\frac{c_1-b_1y^{(k)}}{a_1}
$$

$$
y^{(k+1)}=\frac{c_2-a_2x^{(k)}}{b_2}
$$

Define:

$$
p=\frac{c_1}{a_1},\qquad q=-\frac{b_1}{a_1}
$$

$$
r=\frac{c_2}{b_2},\qquad s=-\frac{a_2}{b_2}
$$

The equations become:

$$
x^{(k+1)}=p+qy^{(k)}
$$

$$
y^{(k+1)}=r+sx^{(k)}
$$

## 3. Derivation of the Scalar Recurrence

From the second equation:

$$
y^{(k)}=r+sx^{(k-1)}
$$

Substitute into the first equation:

$$
x^{(k+1)}=p+q(r+sx^{(k-1)})
$$

Expanding:

$$
x^{(k+1)}=p+qr+qsx^{(k-1)}
$$

Let:

$$
A=p+qr
$$

$$
B=qs
$$

Therefore:

$$
\boxed{x^{(k+1)}=A+Bx^{(k-1)}}
$$

This recurrence is the foundation of the direct-iteration approach.

It connects iterations that are two steps apart, rather than requiring every immediately preceding iteration.

## 4. Derivation of the Direct Iteration Formula

### 4.1 Even-numbered iterations

Starting with:

$$
x^{(2)}=A+Bx^{(0)}
$$

$$
x^{(4)}=A+B(A+Bx^{(0)})
$$

Expanding:

$$
x^{(4)}=A(1+B)+B^2x^{(0)}
$$

Continuing the pattern:

$$
x^{(2m)}=A(1+B+B^2+\cdots+B^{m-1})+B^mx^{(0)}
$$

Using the geometric-series formula:

$$
\boxed{
x^{(2m)}=
A\left(\frac{1-B^m}{1-B}\right)+B^mx^{(0)}
}
$$

### 4.2 Odd-numbered iterations

Similarly:

$$
x^{(1)}=x^{(1)}
$$

$$
x^{(3)}=A+Bx^{(1)}
$$

$$
x^{(5)}=A+B(A+Bx^{(1)})
$$

Therefore:

$$
\boxed{
x^{(2m+1)}=
A\left(\frac{1-B^m}{1-B}\right)+B^mx^{(1)}
}
$$

These expressions apply when $B\neq1$.

When $B=1$, the corresponding expressions are:

$$
\boxed{x^{(2m)}=mA+x^{(0)}}
$$

$$
\boxed{x^{(2m+1)}=mA+x^{(1)}}
$$

The second variable can be recovered using:

$$
\boxed{y^{(n)}=r+sx^{(n-1)}}
$$

for $n\geq1$.

## 5. Worked Example

Consider the Jacobi equations:

$$
x^{(k+1)}=\frac{9-y^{(k)}}{4}
$$

$$
y^{(k+1)}=\frac{8-x^{(k)}}{3}
$$

Initial estimates:

$$
x^{(0)}=0,\qquad y^{(0)}=0
$$

Identify the constants:

$$
p=\frac94,\quad q=-\frac14
$$

$$
r=\frac83,\quad s=-\frac13
$$

Calculate:

$$
A=p+qr
$$

$$
A=\frac94-\frac14\left(\frac83\right)
=\frac{19}{12}
$$

And:

$$
B=qs
$$

$$
B=\left(-\frac14\right)\left(-\frac13\right)
=\frac1{12}
$$

The scalar recurrence becomes:

$$
\boxed{x^{(k+1)}=\frac{19}{12}+\frac1{12}x^{(k-1)}}
$$

The first iteration is:

$$
x^{(1)}=\frac94
$$

### 5.1 Calculate the fifth iteration directly

Since:

$$
5=2(2)+1
$$

Use the odd-iteration formula with $m=2$:

$$
x^{(5)}=
A\left(\frac{1-B^2}{1-B}\right)+B^2x^{(1)}
$$

Substituting:

$$
x^{(5)}
=
\frac{19}{12}
\left(\frac{1-(1/12)^2}{1-1/12}\right)
+ \left(\frac1{12}\right)^2\frac94
$$

Therefore:

$$
\boxed{x^{(5)}=\frac{997}{576}\approx1.730903}
$$

### 5.2 Calculate the corresponding fifth iteration of y

First, calculate $x^{(4)}$ directly:

$$
x^{(4)}=
\frac{19}{12}
\left(\frac{1-(1/12)^2}{1-1/12}\right)
$$

$$
x^{(4)}=\frac{247}{144}
$$

Then:

$$
y^{(5)}=\frac83-\frac13x^{(4)}
$$

$$
\boxed{y^{(5)}=\frac{905}{432}\approx2.094907}
$$

Hence:

$$
\boxed{
X^{(5)}=
\begin{bmatrix}
1.730903\\
2.094907
\end{bmatrix}
}
$$

The fifth iteration is obtained without sequentially calculating iterations two, three and four.

## 6. Algorithm Summary

For a two-variable Jacobi system:

1. Rearrange both equations into Jacobi form.
2. Identify $p$, $q$, $r$, and $s$.
3. Calculate $A=p+qr$ and $B=qs$.

4. Calculate the first iteration.
5. Determine whether the required iteration is even or odd.
6. Apply the corresponding direct formula.
7. Calculate the other variable using $y^{(n)}=r+sx^{(n-1)}$.

## 7. Advantages

- Avoids sequential calculation of every intermediate iteration.
- Reduces repetitive substitution.
- Provides closed-form expressions for selected iterations.
- Connects iterative numerical methods with recurrence relations.
- Can be implemented using a scientific calculator or computer program.

## 8. Limitations

- The derivation presented here is specifically for two-variable Jacobi systems.
- The method does not automatically improve the convergence rate of Jacobi iteration.
- The calculation of the second variable requires the preceding value of the first variable.
- Higher-dimensional systems require more advanced recurrence or matrix-power techniques.
- A direct formula is not necessarily computationally cheaper for every problem.

## 9. Future Development

Possible extensions include:

- Extending the recurrence approach to three-variable Jacobi systems.
- Developing a general matrix-based direct-iteration algorithm.
- Implementing the method in Python.
- Comparing computational efficiency against conventional Jacobi iteration.
- Investigating the relationship between the recurrence coefficient and convergence.

## 10. Conclusion

This project demonstrates that the conventional two-variable Jacobi iteration can be algebraically transformed into a scalar recurrence relation.

The resulting recurrence enables direct computation of selected iterations through a closed-form geometric-series expression.

It provides an alternative way of organizing Jacobi calculations and demonstrates how recurrence relations can simplify repeated numerical computations.

---

**Project status:** Mathematical derivation and worked example.

**Research note:** This project presents an algebraic shortcut derived from the classical Jacobi method. It does not claim that the underlying recurrence relations or geometric-series techniques are new.
