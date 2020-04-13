#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Traduction python des codes Matlab de R. Aster, B. Borchers, C. Thurber

Created on Tue Apr  3 20:03:51 2018

@author: giroux
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl


def kac(A, b, tolx, maxiter, dense=True):
    """
    % Parameter Estimation and Inverse Problems, 2nd edition, 2011
    % by R. Aster, B. Borchers, C. Thurber
    % x=kac(A,b,tolx,maxiter)
    %
    % Implements Kaczmarz's algorithm to solve a system of equations iteratively.
    %
    % Input Parameters:
    %   A - Constraint matrix.
    %   b - right hand side.
    %   tolx - difference tolerence for successive iterations (stopping criteria).
    %   maxiter - maximum iterations (stopping criteria).
    %
    % Output Parameters:
    %      x - solution.
    function x=kac(A,b,tolx,maxiter)

    % First, find the size of the matrix.
    [m,n]=size(A);

    % Make a copy of A' to speed up some accesses.
    AT=A';

    % Setup an initial solution of all zeros.
    x=zeros(n,1);
    iter=0;

    %  Precompute the row norms squared.
    n2=zeros(m,1);
    for i=1:m
      n2(i)=norm(AT(:,i))^2;
    end

    % The main loop performs iterations of Kaczmarz algorithm until
    % maxiters is exceeded or successive iterates differ by less
    % than tolx.
    while (iter <= maxiter)
      % Update the iteration count.
      iter=iter+1;

      %  Start the update cycle with the current solution.
      newx=x;

      %  Perform a cycle of m updates.
      for i=1:m
        newx=newx-((newx'*AT(:,i)-b(i))/(n2(i)))*AT(:,i);
      end

      %  Check for convergence to fixed solution.
      if (norm(newx-x)/(1+norm(x)) < tolx)
        x=newx;
        return;
      end

      %  Update x for the next major iteration.
      x=newx;
    end

    % If no convergence
    disp('Max iterations exceeded.');
    """
    m, n = A.shape
    x = np.zeros((n, 1))
    newx = x.copy()
    n2 = np.zeros((m, 1))
    it = 0
    if dense:
        A = A.toarray()
        AT = A.T
        for i in range(m):
            n2[i] = np.linalg.norm(AT[:, i])**2
        while it <= maxiter:
            it += 1
            for i in range(m):
                newx +=
            if np.linalg.norm(newx-x)/(1+np.linalg.norm(x)) < tolx:
                return newx
            x = newx.copy()
    else:
        AT = A.T
        for i in range(m):
            n2[i] = spl.norm(AT[:, i])**2
        while it <= maxiter:
            it += 1
            for i in range(m):
                newx +=
            if np.linalg.norm(newx-x)/(1+np.linalg.norm(x)) < tolx:
                return newx
            x = newx.copy()

    print('Max iterations exceeded.')
    return x


def art(A, b, tolx, maxiter, dense=True):
    """
    % Parameter Estimation and Inverse Problems, 2nd edition, 2011
    % by R. Aster, B. Borchers, C. Thurber
    % x=art(A,b,tolx,maxiter)
    %
    % Implements the ART algorithm to solve a system of equations iteratively.
    %
    % Input Parameters:
    %   A       - Constraint matrix.
    %   b       - right hand side.
    %   tolx    - difference tolerance for successive iterations (stopping criteria)
    %   maxiter - maximum iterations (stopping criteria).
    %
    % Output Parameters:
    %   x - solution.
    function x=art(A,b,tolx,maxiter)

    % Alpha is a damping factor.  If alpha<1, then we won't take full steps
    % in the ART direction.  Using a smaller value of alpha (say alpha=.75)
    % can help with convergence on some problems.
    alpha=1.0;

    % First, get the size of A.
    [m,n]=size(A);

    % In the A1 array, we convert all nonzero entries in A to 1.
    A1=(A>0);

    % Get transposed copies of the arrays for faster access.
    AP=A';
    A1P=A1';

    %  Precompute N(i) and L(i) factors.
    N=zeros(m,1);
    L=zeros(m,1);
    for i=1:m
      N(i)=sum(A1(i,:));
      L(i)=sum(A(i,:));
    end

    % Start with the zero solution.
    x=zeros(n,1);

    % Start the iteration count at 0.
    iter=0;

    % Now, the main loop.
    while (true)
      % Check to make sure that we haven't exceeded maxiter.
      iter=iter+1;
      if (iter > maxiter)
          disp('Max iterations exceeded.');
          x=newx;
          return;
      end

      % Start the next round of updates with the current solution.
      newx=x;

      % Now, update each of the m constraints.
      for i=1:m
        %  Compute the weighted sum for constraint i.
        q=A1P(:,i)'*newx;

        % We use the more accurate formula for delta.
        delta=b(i)/L(i)-q/N(i);

        % This alternative formula is less accurate and doesn't work nearly as well.
        %
        %    delta=(b(i)-q)/N(i);
        %

        % Now do the update.
        newx=newx+alpha*delta*A1P(:,i);
      end

      % Check for convergence
      if (norm(newx-x)/(1+norm(x)) < tolx)
        x=newx;
        return;
      end

      % No convergence, so setup for the next ART iteration.
      x=newx;
    end
    """
    alpha = 1.0

    m, n = A.shape
    N = np.zeros((m, 1))
    L = np.zeros((m, 1))
    x = np.zeros((n, 1))
    newx = x.copy()
    it = 0
    if dense:
        A = A.toarray()
        A1 = A > 0
        A1P = A1.T
        for i in range(m):
            N[i] = np.sum(A1[i, :])
            L[i] = np.sum(A[i, :])
        while it <= maxiter:
            it += 1
            for i in range(m):
                q = A1P[:, i].T.dot(newx)
                delta =
                newx +=
            if np.linalg.norm(newx-x)/(1+np.linalg.norm(x)) < tolx:
                return newx
            x = newx.copy()
    else:
        A1 = A > 0
        A1P = A1.T
        for i in range(m):
            N[i] = np.sum(A1[i, :])
            L[i] = np.sum(A[i, :])
        while it <= maxiter:
            it += 1
            for i in range(m):
                q = A1P[:, i].T.dot(newx)
                delta =
                newx +=
            if np.linalg.norm(newx-x)/(1+np.linalg.norm(x)) < tolx:
                return newx
            x = newx.copy()

    print('Max iterations exceeded.')
    return x


def sirt(A, b, tolx, maxiter, dense=True):
    """
    % Parameter Estimation and Inverse Problems, 2nd edition, 2011
    % by R. Aster, B. Borchers, C. Thurber
    % x=sirt(A,b,tolx,maxiter)
    %
    % Implements the SIRT algorithm to solve a system of equations iteratively.
    %
    % Input Parameters:
    %   A       - Constraint matrix.
    %   b       - right hand side.
    %   tolx    - difference tolerance for successive iterations (stopping criteria)
    %   maxiter - maximum iterations (stopping criteria).
    %
    % Output Parameters:
    %      x - solution.
    %
    function x=sirt(A,b,tolx,maxiter)

    % First, get the size of A.
    [m,n]=size(A);

    % Alpha is a damping factor.  If alpha<1, then we won't take full steps
    % in the SIRT direction.  Using a smaller value of alpha (say alpha=.75)
    % can help with convergence on some problems.
    alpha=1.0;

    % In the A1 array, we convert all nonzero entries in A to +1.
    A1=(A>0);

    % Get transposed copies of the arrays for faster access.
    AT=A';
    A1T=A1';

    % Start with the zero solution.
    x=zeros(n,1);

    %  Precompute N(i) and L(i) factors.
    N=zeros(m,1);
    L=zeros(m,1);
    NRAYS=zeros(n,1);

    for i=1:m
      N(i)=sum(A1T(:,i));
      L(i)=sum(AT(:,i));
    end

    for i=1:n
      NRAYS(i)=sum(A1(:,i));
    end

    % Start the iteration count at 0.
    iter=0;

    % Now, the main loop, don't loop more than maxiter times
    while (iter<=maxiter)
      iter=iter+1;

      % Start the next round of updates with the current solution.
      newx=x;

      %  Now, compute the updates for all of the rays and all cells, and put
      %  them in a vector called deltax.
      deltax=zeros(n,1);
      for i=1:m,
        %  Compute the approximate travel time for ray i.
        q=A1T(:,i)'*newx;

        % We use the following more accurate formula for delta.
        delta=b(i)/L(i)-q/N(i);

        % This formula is less accurate and doesn't work nearly as well.
        %
        %    delta=(b(i)-q)/N(i);
        %

        %  Perform updates for those cells touched by ray i.
        deltax=deltax+delta*A1T(:,i);
      end

      %  Now, add the average of the updates to the old model.
      %  Note the "./" here.  This means that the averaging is done with
      %  respect to the number of rays that pass through a particular cell!
      newx=newx+alpha*deltax./NRAYS;

      % Check for convergence
      if (norm(newx-x)/(1+norm(x)) < tolx)
        x=newx
        return;
      end

      % No convergence, so setup for the next major iteration.
      x=newx;
    end
    disp('Max iterations exceeded.');
    """
    alpha = 1.0
    m, n = A.shape
    x = np.zeros((n, 1))
    N = np.zeros((m, 1))
    L = np.zeros((m, 1))
    NRAYS = np.zeros((n, 1))
    newx = x.copy()
    it = 0
    if dense:
        A = A.toarray()
        A1 = A > 0
        A1T = A1.T
        for i in range(m):
            N[i] = np.sum(A1[i, :])
            L[i] = np.sum(A[i, :])
        for i in range(n):
            NRAYS[i] = np.sum(A1[:, 1])
        while it <= maxiter:
            it += 1
            deltax = np.zeros((n, 1))
            for i in range(m):
                q = A1T[:, i].T.dot(newx)
                delta =
                deltax +=
            newx += alpha*deltax/NRAYS
            if np.linalg.norm(newx-x)/(1+np.linalg.norm(x)) < tolx:
                return newx
            x = newx.copy()
    else:
        A1 = A > 0
        A1T = A1.T
        for i in range(m):
            N[i] = np.sum(A1[i, :])
            L[i] = np.sum(A[i, :])
        for i in range(n):
            NRAYS[i] = np.sum(A1[:, 1])
        while it <= maxiter:
            it += 1
            deltax = np.zeros((n, 1))
            for i in range(m):
                q = A1T[:, i].T.dot(newx)
                delta =
                deltax +=
            newx += alpha*deltax/NRAYS
            if np.linalg.norm(newx-x)/(1+np.linalg.norm(x)) < tolx:
                return newx
            x = newx.copy()

    print('Max iterations exceeded.')
    return x


if __name__ == '__main__':

    import matplotlib.pyplot as plt
    import time

    # Exemple 6.2
    G = np.zeros((94, 256))
    for i in range(16):
        for j in range(i*16, (i+1)*16):
            G[i, j] = 1.
    for i in range(16):
        for j in range(i, 241+i, 16):
            G[i+16, j] = 1.
    for i in range(16):
        for j in range(i+1):
            G[i+32, i+j*15] = np.sqrt(2.0)
    for i in range(15):
        for j in range(15-i):
            G[i+48, (i+2)*16+j*15-1] = np.sqrt(2.0)
    for i in range(16):
        for j in range(i+1):
            G[i+63, 15-i+17*j] = np.sqrt(2.0)
    for i in range(15):
        for j in range(15-i):
            G[i+79, (i+1)*16+17*j] = np.sqrt(2.0)

    G = sp.csc_matrix(G)

    plt.figure()
    plt.imshow(G.todense())
    plt.colorbar()
    plt.show(block=False)

    mtrue = np.zeros((16, 16))
    mtrue[8, 8] = 1
    mtrue[8, 9] = 1
    mtrue[8, 10] = 1
    mtrue[9, 8] = 1
    mtrue[9, 10] = 1
    mtrue[10, 8] = 1
    mtrue[10, 9] = 1
    mtrue[10, 10] = 1
    mtrue[1, 2] = 1
    mtrue[1, 3] = 1
    mtrue[2, 2] = 1
    mtrue[2, 3] = 1

    mtrue = mtrue.reshape(-1, 1)
    dtrue = G.dot(mtrue)

    d = dtrue + 0.1*np.random.randn(dtrue.size, 1)

    t1 = time.process_time()
    mkac = kac(G, d, 0.0, 200)
    tkac = time.process_time()-t1

    plt.figure()
    plt.subplot(121)
    plt.imshow(mtrue.reshape(16, 16).T)
    plt.colorbar()

    plt.subplot(122)
    plt.imshow(mkac.reshape(16, 16).T)
    plt.colorbar()
    plt.show(block=False)

    t1 = time.process_time()
    mart = art(G, d, 0.0, 200)
    tart = time.process_time()-t1

    plt.figure()
    plt.subplot(121)
    plt.imshow(mtrue.reshape(16, 16).T)
    plt.colorbar()

    plt.subplot(122)
    plt.imshow(mart.reshape(16, 16).T)
    plt.colorbar()
    plt.show(block=False)

    t1 = time.process_time()
    msirt = sirt(G, d, 0.0, 200)
    tsirt = time.process_time()-t1

    plt.figure()
    plt.subplot(121)
    plt.imshow(mtrue.reshape(16, 16).T)
    plt.colorbar()

    plt.subplot(122)
    plt.imshow(msirt.reshape(16, 16).T)
    plt.colorbar()
    plt.show()

    tmp = np.linalg.norm(mtrue)
    print('temps kac: ', tkac)
    print('temps art: ', tart)
    print('temps sirt: ', tsirt)
    print('kac relative error  ', np.linalg.norm(mkac-mtrue)/tmp)
    print('art relative error  ', np.linalg.norm(mart-mtrue)/tmp)
    print('sirt relative error ', np.linalg.norm(msirt-mtrue)/tmp)
