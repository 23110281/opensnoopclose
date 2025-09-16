import os
import sys

# largest prime smaller than 65536.
ADLER_MOD = 65521

def calculate_full_adler32(data: bytes | bytearray) -> int:
    """
    Calculates the full Adler-32 checksum for a block of data.
    
    This implementation follows the specific formulas provided in the paper:
    A = (1 + D1 + D2 + ... + Dn) mod 65521
    B = (n*D1 + (n-1)*D2 + ... + Dn + n) mod 65521
    Adler32(D) = A + B * 65536
    
    Note: This is a direct implementation for demonstrating the paper's logic,
    and is different from the more common rolling implementation of Adler-32
    (e.g., zlib.adler32) where B is a sum of A's.
    
    Args:
        data: A bytes or bytearray object representing the data block.

    Returns:
        The 32-bit Adler checksum as an integer.
    """
    n = len(data)
    
    # Calculate A
    # The initial '1' is added to the sum.
    a = (1 + sum(data)) % ADLER_MOD
    
    # Calculate B
    # The initial 'n' is added to the sum.
    b = n
    for i in range(n):
        # D_i corresponds to data[i-1] in 0-based indexing.
        # The weight for D_i is (n - i + 1).
        # For data[j] (0-indexed), the weight is (n - j).
        b += (n - i) * data[i]
    b %= ADLER_MOD
    
    # Combine A and B into the final checksum
    # B is stored in the upper 16 bits, A in the lower 16 bits.
    return a + (b * 65536)

def calculate_incremental_adler32_single_byte(
    old_checksum: int, 
    n: int, 
    i: int, 
    old_byte: int, 
    new_byte: int
) -> int:
    """
    Incrementally computes the new Adler-32 checksum after a single byte is modified.

    This function implements the following formulas:
    A = C mod 65536
    B = C div 65536
    ∆Di = D′i − Di
    A′ = (A + ∆Di) mod 65521
    B′ = (B + (n + 1 − i) × ∆Di) mod 65521
    Adler32Incr = A′ + B′ × 65536

    Args:
        old_checksum (int): The original checksum of the data block.
        n (int): The total number of bytes in the data block.
        i (int): The 1-based index of the byte that was modified.
        old_byte (int): The original value of the byte (0-255).
        new_byte (int): The new value of the byte (0-255).

    Returns:
        The new 32-bit Adler checksum as an integer.
    """
    if not (1 <= i <= n):
        raise ValueError("Index i must be 1-based and within the bounds of the data length n.")

    # Extract old A and B from the original checksum C
    a = old_checksum % 65536
    b = old_checksum // 65536  # Integer division
    
    # Calculate the difference (delta)
    delta_d = new_byte - old_byte
    
    # Calculate new A'
    new_a = (a + delta_d) % ADLER_MOD
    
    # Calculate new B'
    # The weight (n + 1 - i) corresponds to the 1-based index i.
    weight = n + 1 - i
    new_b = (b + weight * delta_d) % ADLER_MOD
    
    # Recombine into the new checksum
    return new_a + (new_b * 65536)

def calculate_incremental_adler32_contiguous(
    old_checksum: int,
    n: int,
    start_i: int,
    old_bytes: bytes | bytearray,
    new_bytes: bytes | bytearray
) -> int:
    """
    Incrementally recomputes the checksum when multiple contiguous bytes are modified.
    This is an extension of the single-byte modification logic.

    Args:
        old_checksum (int): The original checksum of the data block.
        n (int): The total number of bytes in the data block.
        start_i (int): The 1-based starting index of the modification.
        old_bytes (bytes): The original contiguous bytes being overwritten.
        new_bytes (bytes): The new contiguous bytes.

    Returns:
        The new 32-bit Adler checksum as an integer.
    """
    if len(old_bytes) != len(new_bytes):
        raise ValueError("Old and new byte slices must have the same length.")
    
    # Extract old A and B
    a = old_checksum % 65536
    b = old_checksum // 65536

    # The total change in A is the sum of all individual byte changes.
    delta_a_total = sum(new_bytes) - sum(old_bytes)
    
    delta_b_total = 0
    for j in range(len(old_bytes)):
        # The 1-based index of the current byte being processed
        current_i = start_i + j
        weight = n + 1 - current_i
        delta_d = new_bytes[j] - old_bytes[j]
        delta_b_total += weight * delta_d
        
    new_a = (a + delta_a_total) % ADLER_MOD
    new_b = (b + delta_b_total) % ADLER_MOD
    
    return new_a + (new_b * 65536)


def run_corruption_test(file_path: str):
    """
    Runs a comprehensive series of data modification tests on a given file.
    If the file does not exist, it will be created with sample content.
    """
    try:
        # --- Setup: Ensure the file exists, then read it ---
        if not os.path.exists(file_path):
            print(f"[*] File '{file_path}' not found. Creating it with sample content.")
            # Create a file with enough content for robust testing
            content = ("The quick brown fox jumps over the lazy dog. " * 50 + "\n") * 10
            with open(file_path, "w") as f:
                f.write(content)

        # --- Read the file content into a mutable bytearray ---
        with open(file_path, "rb") as f:
            original_data = bytearray(f.read())
        
        n_initial = len(original_data)
        if n_initial == 0:
            print("[!] Warning: The file is empty. No tests can be performed.")
            return

        print(f"[*] Read {n_initial} bytes from the file '{file_path}' into memory.")
        
        # --- 1. Calculate the initial checksum ---
        initial_checksum = calculate_full_adler32(original_data)
        print(f"\n[Initial State]")
        print(f"  > Full Checksum: {initial_checksum} (0x{initial_checksum:08x})")
        print("=" * 60)

        # --- Test Group A: Size remains constant ---

        # Test A.1: Content Replaced (Contiguous Slice)
        print("[Test A.1: Content Replaced - Size Constant]")
        if n_initial > 50:
            start_index_0based = 50
            start_index_1based = start_index_0based + 1
            length = 16
            
            old_bytes_slice = original_data[start_index_0based : start_index_0based + length]
            new_bytes_slice = os.urandom(length)

            print(f"[*] Replacing {length} bytes at index {start_index_1based}.")
            
            incremental_checksum = calculate_incremental_adler32_contiguous(
                initial_checksum, n_initial, start_index_1based, old_bytes_slice, new_bytes_slice
            )
            
            modified_data = original_data.copy()
            modified_data[start_index_0based : start_index_0based + length] = new_bytes_slice
            full_recalculated_checksum = calculate_full_adler32(modified_data)

            print(f"  > Incremental Checksum: {incremental_checksum:08x}")
            print(f"  > Full Recalculation:   {full_recalculated_checksum:08x}")
            if incremental_checksum == full_recalculated_checksum and full_recalculated_checksum != initial_checksum:
                print("  > \033[92mSUCCESS\033[0m: Incremental update is correct and change was detected.")
            else:
                print("  > \033[91mFAILURE\033[0m: Incremental update incorrect or change missed.")
        else:
            print("  > SKIPPED (file too small).")
        print("-" * 60)

        # Test A.2: Single Byte Corruption (Scribble)
        print("[Test A.2: Single Byte Corruption - Size Constant]")
        mod_index_0based = min(150, n_initial - 1)
        mod_index_1based = mod_index_0based + 1
        old_byte_val = original_data[mod_index_0based]
        new_byte_val = (old_byte_val + 42) % 256
        print(f"[*] Corrupting single byte at index {mod_index_1based}.")
        incremental_checksum = calculate_incremental_adler32_single_byte(
            initial_checksum, n_initial, mod_index_1based, old_byte_val, new_byte_val
        )
        modified_data = original_data.copy()
        modified_data[mod_index_0based] = new_byte_val
        full_recalculated_checksum = calculate_full_adler32(modified_data)
        print(f"  > Incremental Checksum: {incremental_checksum:08x}")
        print(f"  > Full Recalculation:   {full_recalculated_checksum:08x}")
        if incremental_checksum == full_recalculated_checksum and full_recalculated_checksum != initial_checksum:
            print("  > \033[92mSUCCESS\033[0m: Incremental update is correct and change was detected.")
        else:
            print("  > \033[91mFAILURE\033[0m: Incremental update incorrect or change missed.")
        print("-" * 60)

        # Test A.3: Bit-Flip Corruption
        print("[Test A.3: Bit-Flip Corruption - Size Constant]")
        bit_flip_index_0based = min(200, n_initial - 1)
        bit_flip_index_1based = bit_flip_index_0based + 1
        original_byte = original_data[bit_flip_index_0based]
        flipped_byte = original_byte ^ (1 << 3) # Flip the 4th bit
        print(f"[*] Flipping one bit at index {bit_flip_index_1based}.")
        incremental_checksum = calculate_incremental_adler32_single_byte(
            initial_checksum, n_initial, bit_flip_index_1based, original_byte, flipped_byte
        )
        modified_data = original_data.copy()
        modified_data[bit_flip_index_0based] = flipped_byte
        full_recalculated_checksum = calculate_full_adler32(modified_data)
        print(f"  > Incremental Checksum: {incremental_checksum:08x}")
        print(f"  > Full Recalculation:   {full_recalculated_checksum:08x}")
        if incremental_checksum == full_recalculated_checksum and full_recalculated_checksum != initial_checksum:
            print("  > \033[92mSUCCESS\033[0m: Incremental update is correct and change was detected.")
        else:
            print("  > \033[91mFAILURE\033[0m: Incremental update incorrect or change missed.")
        print("-" * 60)

        # --- Test Group B: Size Increases ---

        # Test B.1: New Slice Inserted
        print("[Test B.1: Slice Inserted - Size Increased]")
        if n_initial > 20:
            insert_index = 20
            insert_data = b"***INSERTED DATA***"
            print(f"[*] Inserting {len(insert_data)} bytes at index {insert_index}.")
            print("  > NOTE: Incremental update is NOT possible as file length changes.")
            modified_data = original_data.copy()
            modified_data[insert_index:insert_index] = insert_data
            full_recalculated_checksum = calculate_full_adler32(modified_data)
            print(f"  > Full Recalculation on new data: {full_recalculated_checksum:08x}")
            if full_recalculated_checksum != initial_checksum:
                print("  > \033[92mSUCCESS\033[0m: Change was detected by full recalculation.")
            else:
                print("  > \033[91mFAILURE\033[0m: Change was not detected.")
        else:
            print("  > SKIPPED (file too small).")
        print("-" * 60)

        # Test B.2: New Content Appended
        print("[Test B.2: Content Appended - Size Increased]")
        append_data = b"***THIS IS APPENDED CONTENT***"
        print(f"[*] Appending {len(append_data)} bytes to the end of the file.")
        print("  > NOTE: Incremental update is NOT possible as file length changes.")
        modified_data = original_data.copy()
        modified_data.extend(append_data)
        full_recalculated_checksum = calculate_full_adler32(modified_data)
        print(f"  > Full Recalculation on new data: {full_recalculated_checksum:08x}")
        if full_recalculated_checksum != initial_checksum:
            print("  > \033[92mSUCCESS\033[0m: Change was detected by full recalculation.")
        else:
            print("  > \033[91mFAILURE\033[0m: Change was not detected.")
        print("-" * 60)

        # Test B.3: Replace with Larger Slice
        print("[Test B.3: Replace with Larger Slice - Size Increased]")
        if n_initial > 200:
            replace_start = 180
            replace_len = 10
            new_data = b"This is a much larger slice of data"
            print(f"[*] Replacing {replace_len} bytes with {len(new_data)} bytes at index {replace_start}.")
            print("  > NOTE: Incremental update is NOT possible as file length changes.")
            modified_data = original_data.copy()
            modified_data[replace_start : replace_start + replace_len] = new_data
            full_recalculated_checksum = calculate_full_adler32(modified_data)
            print(f"  > Full Recalculation on new data: {full_recalculated_checksum:08x}")
            if full_recalculated_checksum != initial_checksum:
                print("  > \033[92mSUCCESS\033[0m: Change was detected by full recalculation.")
            else:
                print("  > \033[91mFAILURE\033[0m: Change was not detected.")
        else:
            print("  > SKIPPED (file too small).")
        print("-" * 60)

        # --- Test Group C: Size Decreases ---

        # Test C.1: Slice Removed
        print("[Test C.1: Slice Removed - Size Decreased]")
        if n_initial > 100:
            remove_start = 80
            remove_length = 20
            print(f"[*] Removing {remove_length} bytes starting at index {remove_start}.")
            print("  > NOTE: Incremental update is NOT possible as file length changes.")
            modified_data = original_data.copy()
            del modified_data[remove_start : remove_start + remove_length]
            full_recalculated_checksum = calculate_full_adler32(modified_data)
            print(f"  > Full Recalculation on new data: {full_recalculated_checksum:08x}")
            if full_recalculated_checksum != initial_checksum:
                print("  > \033[92mSUCCESS\033[0m: Change was detected by full recalculation.")
            else:
                print("  > \033[91mFAILURE\033[0m: Change was not detected.")
        else:
            print("  > SKIPPED (file too small).")
        print("-" * 60)
        
        # Test C.2: Content Truncated
        print("[Test C.2: File Truncated - Size Decreased]")
        if n_initial > 100:
            truncate_size = n_initial // 2
            print(f"[*] Truncating file to its first {truncate_size} bytes.")
            print("  > NOTE: Incremental update is NOT possible as file length changes.")
            modified_data = original_data[:truncate_size]
            full_recalculated_checksum = calculate_full_adler32(modified_data)
            print(f"  > Full Recalculation on new data: {full_recalculated_checksum:08x}")
            if full_recalculated_checksum != initial_checksum:
                print("  > \033[92mSUCCESS\033[0m: Change was detected by full recalculation.")
            else:
                print("  > \033[91mFAILURE\033[0m: Change was not detected.")
        else:
            print("  > SKIPPED (file too small).")
        print("-" * 60)

        # Test C.3: Replace with Smaller Slice
        print("[Test C.3: Replace with Smaller Slice - Size Decreased]")
        if n_initial > 300:
            replace_start = 250
            replace_len = 30
            new_data = b"smaller"
            print(f"[*] Replacing {replace_len} bytes with {len(new_data)} bytes at index {replace_start}.")
            print("  > NOTE: Incremental update is NOT possible as file length changes.")
            modified_data = original_data.copy()
            modified_data[replace_start : replace_start + replace_len] = new_data
            full_recalculated_checksum = calculate_full_adler32(modified_data)
            print(f"  > Full Recalculation on new data: {full_recalculated_checksum:08x}")
            if full_recalculated_checksum != initial_checksum:
                print("  > \033[92mSUCCESS\033[0m: Change was detected by full recalculation.")
            else:
                print("  > \033[91mFAILURE\033[0m: Change was not detected.")
        else:
            print("  > SKIPPED (file too small).")
        print("=" * 60)

    except FileNotFoundError:
        print(f"Error: Could not find the file at {file_path}", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
    finally:
        # The script no longer deletes the file, so the user can inspect it.
        print(f"\n[*] Test run finished for '{file_path}'.")


if __name__ == "__main__":
    # You can change this to any file you want to test.
    # If it doesn't exist, the script will create it.
    FILE_PATH = "test.txt"

    run_corruption_test(FILE_PATH)
