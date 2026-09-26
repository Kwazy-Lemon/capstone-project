# Python-Based SPICE-Like Circuit Simulator

## Group 24

McGill University — Capstone Project

---

## Project Title

**Python-Based SPICE-Like Circuit Simulator**

---

## Group Members

| # | Name | Student ID | Email |
|---|---|---|---|
| 1 | Jianhao Wu | 261073381 | jianhao.wu@mail.mcgill.ca |
| 2 | Zhuoyao Luo | 261139787 | zhuoyao.luo@mail.mcgill.ca |
| 3 | Minho Chang | 261163554 | minho.chang@mail.mcgill.ca |
| 4 | Jiahao Zhang | 261138548 | jiahao.zhang3@mail.mcgill.ca |
| 5 | Zhihang Chen | 261070033 | zhihang.chen@mail.mcgill.ca |

---

## Project Advisor

**Prof. Roni Khazaka**

Associate Professor  
Associate Dean (Academic Programs)  
Department of Electrical and Computer Engineering  
McGill University

**Email:** roni.khazaka@mcgill.ca

**Research Area:** Integrated Circuits and Systems

---

## Intellectual Property

The intellectual property and source code resulting from this project will be jointly owned by the faculty advisor and all student team members, subject to applicable McGill University policies and any applicable project agreements.

All project participants will retain the ability to further develop and use the resulting work for appropriate academic, research, educational, and other permitted purposes.

---

## Non-Disclosure Agreement (NDA)

No NDA is currently required for this project.

The project is an academic software development project and does not currently involve confidential or proprietary industrial information.

If confidential or proprietary resources are introduced during the project, the team will follow the applicable NDA and confidentiality requirements.

---

## Group Meetings and Meetings with Advisor

The five team members will meet regularly to coordinate project development, review progress, discuss technical issues, and integrate individual contributions.

The team will also meet regularly with the project advisor, Prof. Roni Khazaka, to:

- Present project progress
- Discuss technical questions
- Receive feedback
- Review major project decisions

Meetings will primarily be conducted online through Microsoft Teams.

The regular meeting window is:

**Tuesday, 3:30 p.m. – 5:00 p.m.**

Approximately 30 minutes will be allocated for each meeting.

The team has established a GitHub repository for collaborative development.

GitHub will be used for:

- Source-code management
- Version control
- Task tracking
- Issue management
- Code review
- Documentation
- Team coordination

---

## Project Requirements

The project is primarily a software-based circuit simulation project.

The final simulator will be developed in **Python** and will use a **Python API** for programmatic circuit construction.

The project will require:

- Python
- NumPy
- SciPy
- Matplotlib
- Git and GitHub
- Standard computers and a Python development environment
- Numerical libraries supporting matrix and sparse-matrix computation

No specialized physical laboratory equipment or hardware is currently required.

### MATLAB

MATLAB will also be used during the early stage of the project for preliminary numerical experiments and validation of MNA formulations.

MATLAB will serve as a prototyping and learning tool rather than the implementation platform for the final simulator.

---

## Project Abstract

This project aims to develop a Python-based SPICE-like circuit simulator that allows users to construct and simulate electrical circuits programmatically through a Python API.

Modified Nodal Analysis (MNA) will be used as the fundamental mathematical formulation for constructing circuit equations and matrices.

The simulator will support:

- DC analysis
- AC/frequency-domain analysis
- Transient/time-domain analysis

The project is motivated by the flexibility of an API-based approach.

Instead of relying on a traditional SPICE netlist, users will be able to construct circuits directly within Python programs, making the simulator suitable for:

- Parameter sweeps
- Optimization
- Automated design exploration
- Simulation-data generation
- Future integration with other computational tools

The project will initially focus on linear circuit elements and systematic MNA matrix construction using established numerical libraries and sparse-matrix techniques.

The framework will then be extended to selected nonlinear devices, including:

- Diodes
- CMOS/MOSFETs
- BJTs

Expected outcomes include:

- A functional Python-based circuit simulation framework
- MNA and matrix-assembly modules
- DC/AC/transient analysis capabilities
- Selected device models
- Validation tests
- Technical documentation

---

## Project Status

This repository contains the source code and documentation for **Group 24's Python-Based SPICE-Like Circuit Simulator**.

Development will proceed as the project progresses.

---
