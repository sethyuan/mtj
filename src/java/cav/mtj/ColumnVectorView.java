package cav.mtj;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.AbstractVector;
import no.uib.cipr.matrix.DenseVector;

/**
 * Vector backed by a matrix.
 */
public class ColumnVectorView extends AbstractVector {
  private Matrix m;
  private int column;

  public ColumnVectorView(Matrix m, int column) {
    super(m.numRows());
    this.m = m;
    this.column = column;
  }

  @Override
  public void add(int row, double value) {
    m.add(row, column, value);
  }

  @Override
  public DenseVector copy() {
    return new DenseVector(this);
  }

  @Override
  public double get(int row) {
    return m.get(row, column);
  }

  @Override
  public void set(int row, double value) {
    m.set(row, column, value);
  }
}
