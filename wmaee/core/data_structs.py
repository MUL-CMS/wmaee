class DotDict(dict):
    """
    A custom dictionary class that allows accessing items using dot notation.

    Example:
    ```
    my_dict = DotDict({"item": "value"})
    item_value = my_dict.item
    print(item_value)  # Output: "value"
    ```
    """

    def __getattr__(self, attr):
        """
        Retrieve a dictionary item using dot notation.

        Args:
            attr (str): The key for the item to retrieve.

        Returns:
            Any: The value associated with the given key.

        Raises:
            AttributeError: If the key does not exist in the dictionary.

        Example:
        ```
        my_dict = DotDict({"item": "value"})
        item_value = my_dict.item
        ```
        """
        if attr in self:
            return self[attr]
        raise AttributeError(f"'DotDict' object has no attribute '{attr}'")