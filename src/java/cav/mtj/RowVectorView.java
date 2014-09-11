package cav.mtj;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.AbstractVector;
import no.uib.cipr.matrix.DenseVector;

/**
 * Vector backed by a matrix.
 */
public class RowVectorView extends AbstractVector {
  private Matrix m;
  private int row;

  public RowVectorView(Matrix m, int row) {
    super(m.numColumns());
    this.m = m;
    this.row = row;
  }

  @Override
  public void add(int column, double value) {
    m.add(row, column, value);
  }

  @Override
  public DenseVector copy() {
    return new DenseVector(this);
  }

  @Override
  public double get(int column) {
    return m.get(row, column);
  }

  @Override
  public void set(int column, double value) {
    m.set(row, column, value);
  }
}
