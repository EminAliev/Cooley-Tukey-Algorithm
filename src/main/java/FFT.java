

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class FFT {

    private static final Complex ZERO = new Complex(0, 0);

    private FFT() {
    }


    public static Complex[] fft(Complex[] x) {
        int n = x.length;


        if (n == 1) {
            return new Complex[]{x[0]};
        }


        if (n % 2 != 0) {
            throw new IllegalArgumentException("n is not a power of 2");
        }


        Complex[] even = new Complex[n / 2];
        for (int k = 0; k < n / 2; k++) {
            even[k] = x[2 * k];
        }
        Complex[] q = fft(even);


        Complex[] odd = even;  // reuse the array
        for (int k = 0; k < n / 2; k++) {
            odd[k] = x[2 * k + 1];
        }
        Complex[] r = fft(odd);


        Complex[] y = new Complex[n];
        for (int k = 0; k < n / 2; k++) {
            double kth = -2 * k * Math.PI / n;
            Complex wk = new Complex(Math.cos(kth), Math.sin(kth));
            y[k] = q[k].plus(wk.times(r[k]));
            y[k + n / 2] = q[k].minus(wk.times(r[k]));
        }
        return y;
    }


    public static Complex[] ifft(Complex[] x) {
        int n = x.length;
        Complex[] y = new Complex[n];


        for (int i = 0; i < n; i++) {
            y[i] = x[i].conjugate();
        }


        y = fft(y);


        for (int i = 0; i < n; i++) {
            y[i] = y[i].conjugate();
        }


        for (int i = 0; i < n; i++) {
            y[i] = y[i].scale(1.0 / n);
        }

        return y;

    }


    public static Complex[] cconvolve(Complex[] x, Complex[] y) {


        if (x.length != y.length) {
            throw new IllegalArgumentException("Dimensions don't agree");
        }

        int n = x.length;


        Complex[] a = fft(x);
        Complex[] b = fft(y);


        Complex[] c = new Complex[n];
        for (int i = 0; i < n; i++) {
            c[i] = a[i].times(b[i]);
        }


        return ifft(c);
    }

    public static Complex[] convolve(Complex[] x, Complex[] y) {
        Complex[] a = new Complex[2 * x.length];
        for (int i = 0; i < x.length; i++)
            a[i] = x[i];
        for (int i = x.length; i < 2 * x.length; i++)
            a[i] = ZERO;

        Complex[] b = new Complex[2 * y.length];
        for (int i = 0; i < y.length; i++)
            b[i] = y[i];
        for (int i = y.length; i < 2 * y.length; i++)
            b[i] = ZERO;

        return cconvolve(a, b);
    }


    private static void show(Complex[] x, String title) {
        System.out.println(title);
        System.out.println("-------------------");
        for (int i = 0; i < x.length; i++) {
            System.out.println(x[i]);
        }
        System.out.println();
    }


    public static void main(String[] args) throws FileNotFoundException {
        Scanner scanner = new Scanner(new File("001.dat"));
        int[] tall = new int[256];
        int j = 0;
        while (scanner.hasNextInt()) {
            tall[j++] = scanner.nextInt();
        }
        Complex[] cinput = new Complex[tall.length];
        for (int i = 0; i < tall.length; i++)
            cinput[i] = new Complex(tall[i], 0.0);


        Complex[] y = fft(cinput);
        show(y, "y = fft(x)");


        Complex[] z = ifft(cinput);
        show(z, "z = ifft(y)");

        Complex[] c = cconvolve(cinput, cinput);
        show(c, "c = cconvolve(x, x)");

        Complex[] d = convolve(cinput, cinput);
        show(d, "d = convolve(x, x)");
    }

}
