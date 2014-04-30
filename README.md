# Sparse and group regression models in Portfolio Optimization 
## Introduction

This repo contains the implementation of models studied, analysed and proposed in ["The effects of Sparse and Group Regression models in Portfolio Optimization"](../blob/master/paper.pdf?raw=true).

This implementation focuses on finding the effects of Sparse and Group regression approaches to portfolio optimization problems in finance.

## Motivation

Current approaches to portfolio optimization consider stocks as individual entities, and do not exploit the grouping/classifying information available (e.g. Financial Sectors, Industries, Type, etc).

This paper proposes a novel approach to Index Tracking - namely, a sparse, group and sparse group approach.

## Implementation
This repo contains the implementation of the following models:

#### Feature Regression Models
* Absolute Values
* Conditional-Value-at-Risk (CVaR) Optimization
* Norm-Constrained CVaR Optimization
* Lasso

#### Group Regression Models
* Group Selection
* Group Lasso
* Sparse Group Lasso

## Requirements
This implementation requires the [CVX library](http://cvxr.com/cvx/download/) for solving the convex optimization problems.

## Usage
Tests were built to provide intuition when implementing the Sparse Group Regression model into a set of data, however, understanding on these models is required for an effective use;
