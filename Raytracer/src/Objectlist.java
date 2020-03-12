/**
 * Klasse implementiert eine Liste von Objekten
 */
public class Objectlist {
	public int N; // Anzahl der Objekte in der Liste
	public Listelement first;
	public Listelement last;
	private Listelement actual;

	public Objectlist() {
		N = 0;
		first = null;
		actual = first;
		last = null;
	}

	public void rewind() {
		actual = first;
	}

	public void next() {
		actual = actual.next;
	}

	public Object getActual() {
		if (actual != null) {
			return actual.data;
		}
		else {
			return null;
		}
	}

	public void add(Object obj) {
		final Listelement e = new Listelement(obj);
		if (N == 0) {
			first = e;
			last = first;

		}
		else {
			last.next = e;
			e.prev = last;
			last = e;
		}
		N++;
	}

	public void add(Objectlist list) {
		if (list.N > 0) {
			if (N == 0) {
				this.first = list.first;
				this.last = list.last;
			}
			else {
				this.last.next = list.first;
				list.first.prev = this.last;
				this.last = list.last;
			}
			N += list.N;
		}
	}

	public Object[] getArray() {
		Object[] array;
		if (N > 0) {
			array = new Object[N];
			actual = first;
			for (int i = 0; i < N; i++) {
				array[i] = actual.data;
				actual = actual.next;
			}
		}
		else {
			array = null;
		}
		return array;
	}

	class Listelement {
		Object data;
		Listelement next;
		Listelement prev;

		Listelement(Object d) {
			data = d;
			next = null;
			prev = null;
		}
	}
}
