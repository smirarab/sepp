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

package edu.rice.cs.bioinfo.programs.phylonet.structs;

import java.util.Arrays;
import java.util.BitSet;
import java.util.Iterator;



/**
 * This class provides an implementation of an arbitrarily long bit vector.  It is designed
 * to minimize space usage.
 * 
 * @author Derek Ruths
 *
 */
public class BitVector implements Iterable<Boolean> {

	// constants
	public static final int BITSTRING_FORMAT = 0;
	
	private static final int BLOCK_SIZE = Integer.SIZE;
	
	// fields
	/**
	 * This field specifies how many bits are encoded in this vector.
	 */
	protected int _length;
	
	/**
	 * Each int holds {@link BLOCK_SIZE} bits.  In order to allow for the vector
	 * to hold more than <code>BLOCK_SIZE</code> bits, multiple ints must be used.
	 */
	protected int[] _blocks;
	protected int _last_block_mask;
	
	// constructors
	/**
	 * Creates a zero-ed out bit vector of length <code>len</code>.
	 * 
	 * @param len is the number of bits this vector holds.
	 */
	public BitVector(int len) {
		_length = len;
		int num_ints = getBlock(len - 1) + 1;
		
		_blocks = new int[num_ints];
		
		// setup the last block mask
		_last_block_mask = 0x0000;
		
		int last_pos = _length - (_blocks.length - 1) * BLOCK_SIZE;
		
		for(int i = 0; i < last_pos; i++) {
			_last_block_mask |= (0x0001 << i);
		}
	}
	
	/**
	 * Creates a bit vector with the same length and value as the specified bit vector.
	 */
	public BitVector(BitVector bv) {
		this(bv.getLength());
		
		setValue(bv);
	}
	
	// methods
	/**
	 * This method returns which integer in the <code>_blocks</code> array the specified
	 * position is found in.
	 */
	protected int getBlock(int pos) {
		return (int) Math.floor((double) pos / (double) BLOCK_SIZE);
	}
	
	/**
	 * Increment the value of this bit-vector
	 */
	public void increment() {
		boolean carry = true;
		int pos = 0;
		
		while(carry && pos < _length) {
			
			if(getValue(pos)) {
				setValue(pos,false);
				carry = true;
			} else {
				setValue(pos,true);
				carry = false;
			}
			
			pos++;
		}
		
		return;
	}

	/**
	 * Decrement the value of this bit-vector
	 */
	public void decrement() {
		
		int pos = 0;
		
		// find the least-significant position with value = true
		while(pos < _length && !getValue(pos)) { pos++; }
		
		// subtract that and spread the ones out
		if(pos < _length) {
			setValue(pos, false);
		}
		
		pos--;
		
		while(pos >= 0) {
			setValue(pos--, true);
		}
		
		return;
	}
	
	/**
	 * @return the value of a single bit in the vector.
	 */
	public boolean getValue(int pos) {
		int block_num = getBlock(pos);
		int offset = pos - (block_num * BLOCK_SIZE);
		
		return ((0x0001 & (_blocks[block_num] >> offset)) != 0x0000);
	}
	
	/**
	 * Set the value of a single bit in this vector.
	 * 
	 * @param pos is the position of the bit to change
	 * @param value is the new value of the bit
	 */
	public void setValue(int pos, boolean value) {
		int block_num = getBlock(pos);
		int offset = pos - (block_num * BLOCK_SIZE);
		
		if(value) {
			_blocks[block_num] |= (0x0001 << offset);
		} else {
			_blocks[block_num] &= (~(0x0001 << offset));
		}
	}
	
	/**
	 * OR the value of this bit vector with another bit vector.
	 * 
	 * @param bv must be the same length as this bit vector.
	 */
	public void or(BitVector bv) {
		if(_length != bv._length) {
			throw new RuntimeException("Unequal bitvector lengths");
		}
		
		for(int i = 0; i < _blocks.length; i++) {
			_blocks[i] |= bv._blocks[i];
		}
		
		// wipe the bits that aren't being used in the last block
		_blocks[_blocks.length-1] &= _last_block_mask;
	}
	
	/**
	 * AND the value of this bit vector with another bit vector.
	 * 
	 * @param bv must be the same length as this bit vector.
	 */
	public void and(BitVector bv) {
		if(_length != bv._length) {
			throw new RuntimeException("Unequal bitvector lengths");
		}
		
		for(int i = 0; i < _blocks.length; i++) {
			_blocks[i] &= bv._blocks[i];
		}
		
		// wipe the bits that aren't being used in the last block
		_blocks[_blocks.length-1] &= _last_block_mask;
	}
	
	/**
	 * XOR the value of this bit vector with another bit vector.
	 * 
	 * @param bv must be the same length as this bit vector.
	 */
	public void xor(BitVector bv) {
		if(_length != bv._length) {
			throw new RuntimeException("Unequal bitvector lengths");
		}
		
		for(int i = 0; i < _blocks.length; i++) {
			_blocks[i] ^= bv._blocks[i];
		}
		
		// wipe the bits that aren't being used in the last block
		_blocks[_blocks.length-1] &= _last_block_mask;
	}
	
	/**
	 * Flip the bits of each bit in this vector
	 */
	public void not() {
		for(int i = 0; i < _blocks.length; i++) {
			_blocks[i] = ~_blocks[i];
		}
		
		// wipe the bits that aren't being used in the last block
		_blocks[_blocks.length-1] &= _last_block_mask;
	}
	
	/**
	 * Set the value of this bit vector to exactly the value of 
	 * the provided vector.
	 * 
	 * @param bv must be a bit vector of the same length.
	 */
	public void setValue(BitVector bv) {
		if(_length != bv._length) {
			throw new RuntimeException("Unequal bitvector lengths");
		}
		
		for(int i = 0; i < _blocks.length; i++) {
			_blocks[i] = bv._blocks[i];
		}
	}
	
	/**
	 * @return the number of bits in this vector
	 */
	public int getLength() {
		return _length;
	}

	/**
	 * @return an iterator that will iterate through the bits of this bit-vector in order of
	 * least-significant to most-significant.
	 */
	public Iterator<Boolean> iterator() {
		return new Iterator<Boolean>() {

			int pos = 0;
			
			public boolean hasNext() {
				return (pos < _length);
			}

			public Boolean next() {
				return getValue(pos++);
			}

			public void remove() {
				throw new RuntimeException("Cannot remove elements from a bit-vector");
			}
		};
	}

	public String toString() {
		char[] val = new char[_length];
		
		for(int i = 0; i < val.length; i++) {
			val[i] = (getValue(i))?'1':'0';
		}
		
		return String.valueOf(val);
	}
	
	public String toString(int format) {
		switch(format) {
		case BITSTRING_FORMAT:
			return toString();
		default:
			throw new RuntimeException("Unknown format " + format);
		}
	}
	
	public int hashCode() {
		// return the hash
		return Arrays.hashCode(_blocks);
	}
	
	/**
	 * @return the number of 1's in this bitvector.
	 */
	public int countOnes() {
		int num = 0;
		
		for(int block : _blocks) {
			num += Integer.bitCount(block);
		}
		
		return num;
	}
	
	public boolean equals(Object obj) {
		if(obj instanceof BitVector) {
			BitVector bv = (BitVector) obj;
			
			if(bv._length != _length) {
				return false;
			}
			
			return Arrays.equals(_blocks,bv._blocks);
		} else {
			return false;
		}
	}
	
	public BitSet toBitSet(){
		BitSet bs = new BitSet(getLength());
		for(int index=0; index<this.getLength(); index++){
			if(this.getValue(index)){
				bs.set(index);
			}
		}
		return bs;
	}
	
	public boolean containsBitVector(BitVector bv){
		if(this.getLength()!=bv.getLength()){
			return false;
		}
		BitVector temp = new BitVector(bv);
		temp.and(this);
		return temp.equals(bv);
	}
}
