# Week 1 — Project Planning and Initial Algorithm Study

## 1. Project Objectives

The main objectives for Week 1 are to further clarify the technical direction of the SPICE-like circuit simulator and establish a foundation for the subsequent Python implementation.

The work will focus on two main areas:

- Developing a **MATLAB-based prototype** to study MNA and basic circuit analysis methods.
- Developing a preliminary **Python software architecture** for the final simulator.

---

## 2. MATLAB-Based Circuit Analysis

A simple circuit will be manually formulated using MNA in MATLAB to study:

- **DC Analysis**
- **AC / Frequency-Domain Analysis**
- **Transient / Time-Domain Analysis**

The prototype will focus on understanding how MNA equations are constructed and solved for different types of circuit analysis.

In parallel, a simple **diode circuit** will be used to investigate nonlinear DC analysis, including:

- Nonlinear function \(F(x)\)
- Jacobian
- Newton Iteration
- Convergence

The purpose of the prototype is to understand the underlying algorithms rather than to develop a complete automated simulator.

---

## 3. Preliminary Python Architecture

Based on the understanding gained from the MATLAB prototype, the team will develop a preliminary architecture for the final Python simulator.

The design will consider:

- Classes and Objects
- Properties and Methods
- Relationships between objects
- Circuit and circuit-element organization
- MNA assembly
- Numerical solving
- DC, AC, and Transient Analysis

The architecture should also consider how both linear and nonlinear circuit elements can be incorporated into a common framework.

---

## 4. Key Technical Considerations

### Node Management

Develop an approach for mapping user-defined node names, such as `Vin`, `out`, and `N1`, to internal integer indices required by MNA.

Python dictionaries/hash tables will be investigated as a possible approach for efficient node management.

### Sparse Matrix Construction

Investigate how sparse matrices can be incorporated from the beginning of the project, including efficient collection and construction of matrix entries for MNA assembly.

---

## 5. Expected Deliverables

By the end of Week 1, the team aims to complete:

1. A simple **MATLAB MNA prototype** for a linear circuit.
2. Initial implementations of **DC, AC, and Transient Analysis**.
3. An initial study of **nonlinear DC analysis using a diode**.
4. A preliminary **Python class and software architecture design**.
5. A list of technical questions and issues to be discussed in the next meeting.
