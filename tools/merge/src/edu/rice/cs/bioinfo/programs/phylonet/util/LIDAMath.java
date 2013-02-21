/*
 * Copyright (c) 2012 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phylonet.util;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.Comparator;

/**
 * This class provides various basic mathematical functions. 
 * 
 * @author Derek Ruths
 */
public class LIDAMath {

	/**
	 * Constructor makes it impossible to instantiate this class
	 */
	private LIDAMath() {}
	
	// static methods
	/**
	 * Computes log<sub>10</sub>(x)
	 */
	public static double log10(double x) {
		return Math.log(x) / Math.log(10);
	}
	
	/**
	 * Computes log<sub>base</sub>(x).
	 * 
	 * @param base is the base of the logarithm
	 * @param x is the value the log is evaluated at
	 * 
	 * @return the value of log<sub>base</sub>(x)
	 */
	public static double log(double base, double x) {
		return Math.log(x) / Math.log(base);
	}
	
	/**
	 * A method to compare two numbers.
	 * 
	 * @param <N> is the type of number being compared
	 */
	public static <N extends Number> int compare(N n1, N n2) {

		if(n1 instanceof Double) {
			Double d1 = (Double) n1;
			Double d2 = (Double) n2;
			
			return d1.compareTo(d2);
		} else if(n1 instanceof Float) {
			Float d1 = (Float) n1;
			Float d2 = (Float) n2;
			
			return d1.compareTo(d2);
		} else if(n1 instanceof Long) {
			Long d1 = (Long) n1;
			Long d2 = (Long) n2;
			
			return d1.compareTo(d2);
		} else if(n1 instanceof Short) {
			Short d1 = (Short) n1;
			Short d2 = (Short) n2;
			
			return d1.compareTo(d2);
		} else if(n1 instanceof Integer) {
			Integer d1 = (Integer) n1;
			Integer d2 = (Integer) n2;
			
			return d1.compareTo(d2);
		} else if(n1 instanceof Byte) {
			Byte d1 = (Byte) n1;
			Byte d2 = (Byte) n2;
			
			return d1.compareTo(d2);
		} else if(n1 instanceof BigInteger) {
			BigInteger d1 = (BigInteger) n1;
			BigInteger d2 = (BigInteger) n2;
			
			return d1.compareTo(d2);
		} else if(n1 instanceof BigDecimal) {
			BigDecimal d1 = (BigDecimal) n1;
			BigDecimal d2 = (BigDecimal) n2;
			
			return d1.compareTo(d2);
		} else {
			throw new RuntimeException("Number type " + n1.getClass() + " is unsupported by MaxHeap");
		}
	}

	public static <N extends Number> Comparator<N> getComparator() {
		return new Comparator<N>() {
			public int compare(N n1, N n2) {
				return LIDAMath.compare(n1,n2);
			}
		};
	}
}
