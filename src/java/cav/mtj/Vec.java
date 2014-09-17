package cav.mtj;

import java.lang.Iterable;
import java.util.Iterator;

import no.uib.cipr.matrix.Vector;

public final class Vec implements Iterable<Double> {
  private class VecIterator implements Iterator<Double> {
    private int size;
    private int index = 0;

    public VecIterator() {
      size = mtj.size();
    }

    @Override
    public boolean hasNext() {
      return index < size;
    }

    @Override
    public Double next() {
      return mtj.get(index++);
    }

    @Override
    public void remove() {
      throw new UnsupportedOperationException();
    }
  }

  Vector mtj;

  public Vec(Vector vector) {
    mtj = vector;
  }

  @Override
  public Iterator<Double> iterator() {
    return new VecIterator();
  }

  public Vec add(double alpha, Vec y) {
    mtj.add(alpha, y.mtj);
    return this;
  }

  public void add(int index, double value) {
    mtj.add(index, value);
  }

  public Vec add(Vec y) {
    mtj.add(y.mtj);
    return this;
  }

  public Vec copy() {
    return new Vec(mtj.copy());
  }

  public double dot(Vec y) {
    return mtj.dot(y.mtj);
  }

  public double get(int index) {
    return mtj.get(index);
  }

  public double norm(Vector.Norm type) {
    return mtj.norm(type);
  }

  public Vec scale(double alpha) {
    mtj.scale(alpha);
    return this;
  }

  public Vec set(double alpha, Vec y) {
    mtj.set(alpha, y.mtj);
    return this;
  }

  public void set(int index, double value) {
    mtj.set(index, value);
  }

  public Vec set(Vec y) {
    mtj.set(y.mtj);
    return this;
  }

  public int size() {
    return mtj.size();
  }

  public Vec zero() {
    mtj.zero();
    return this;
  }

  public String toString() {
    return mtj.toString();
  }

  public Vector mtj() {
    return mtj;
  }
}
