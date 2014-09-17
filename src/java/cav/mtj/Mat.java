package cav.mtj;

import java.lang.Iterable;
import java.util.Iterator;

import no.uib.cipr.matrix.Matrix;

public final class Mat implements Iterable<Double> {
  private class MatIterator implements Iterator<Double> {
    private int numRows, numCols;
    private int row = 0;
    private int col = 0;

    public MatIterator() {
      numRows = mtj.numRows();
      numCols = mtj.numColumns();
    }

    @Override
    public boolean hasNext() {
      return row < numRows && col < numCols;
    }

    @Override
    public Double next() {
      double val = mtj.get(row, col);

      if (col < numCols - 1) {
        col++;
      } else {
        row++;
        col = 0;
      }

      return val;
    }

    @Override
    public void remove() {
      throw new UnsupportedOperationException();
    }
  }

  Matrix mtj;

  public Mat(Matrix matrix) {
    mtj = matrix;
  }

  @Override
  public Iterator<Double> iterator() {
    return new MatIterator();
  }

  public Mat add(double alpha, Mat B) {
    mtj.add(alpha, B.mtj);
    return this;
  }

  public void add(int row, int column, double value) {
    mtj.add(row, column, value);
  }

  public Mat add(Mat B) {
    mtj.add(B.mtj);
    return this;
  }

  public Mat copy() {
    return new Mat(mtj.copy());
  }

  public double get(int row, int column) {
    return mtj.get(row, column);
  }

  public boolean isSquare() {
    return mtj.isSquare();
  }

  public Mat mult(double alpha, Mat B, Mat C) {
    mtj.mult(alpha, B.mtj, C.mtj);
    return C;
  }

  public Vec mult(double alpha, Vec x, Vec y) {
    mtj.mult(alpha, x.mtj, y.mtj);
    return y;
  }

  public Mat mult(Mat B, Mat C) {
    mtj.mult(B.mtj, C.mtj);
    return C;
  }

  public Vec mult(Vec x, Vec y) {
    mtj.mult(x.mtj, y.mtj);
    return y;
  }

  public Mat multAdd(double alpha, Mat B, Mat C) {
    mtj.multAdd(alpha, B.mtj, C.mtj);
    return C;
  }

  public Vec multAdd(double alpha, Vec B, Vec C) {
    mtj.multAdd(alpha, B.mtj, C.mtj);
    return C;
  }

  public Mat multAdd(Mat B, Mat C) {
    mtj.multAdd(B.mtj, C.mtj);
    return C;
  }

  public Vec multAdd(Vec B, Vec C) {
    mtj.multAdd(B.mtj, C.mtj);
    return C;
  }

  public double norm(Matrix.Norm type) {
    return mtj.norm(type);
  }

  public int numColumns() {
    return mtj.numColumns();
  }

  public int numRows() {
    return mtj.numRows();
  }

  public Mat scale(double alpha) {
    mtj.scale(alpha);
    return this;
  }

  public Mat set(double alpha, Mat B) {
    mtj.set(alpha, B.mtj);
    return this;
  }

  public void set(int row, int column, double value) {
    mtj.set(row, column, value);
  }

  public Mat set(Mat B) {
    mtj.set(B.mtj);
    return this;
  }

  public Mat solve(Mat B, Mat X) {
    mtj.solve(B.mtj, X.mtj);
    return X;
  }

  public Vec solve(Vec B, Vec X) {
    mtj.solve(B.mtj, X.mtj);
    return X;
  }

  public Mat transpose() {
    mtj.transpose();
    return this;
  }

  public Mat transpose(Mat B) {
    mtj.transpose(B.mtj);
    return B;
  }

  public Mat zero() {
    mtj.zero();
    return this;
  }

  @Override
  public String toString() {
    return mtj.toString();
  }

  public Matrix mtj() {
    return mtj;
  }
}
