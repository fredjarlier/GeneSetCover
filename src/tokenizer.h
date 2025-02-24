/**
 * \file tokenizer.h
 * \date 3 juil. 2014
 * \brief String tokenizer by char
 * \author Nicolas Fedy
 */

/**
 * \brief Returns the part of the string before the delimiter was found
 * \pre The memory has to be allocated before
 *
 * \param[in] str the string that has to be tokenized
 * \param[in] delim the delimiter to tokenize the string by.
 * \param[out] token first token
 *
 * \return whether 0 or 1 depending on how it went
 */
int tokenizer(char* str, const char delim, char* token);
