/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
  ============================================================================
  Copyright (C) 2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  ============================================================================
  Author: Rene Maerker <rene.maerker@fu-berlin.de>
  06.01.2011
  ============================================================================
  todo enter description of this file
  ==========================================================================*/

#ifndef SCORE_JOURNAL_H_
#define SCORE_JOURNAL_H_

namespace SEQAN_NAMESPACE_MAIN
{
	struct ExtendedScore_;
	typedef Tag<ExtendedScore_> const ExtendedScore;

	template <typename TScoreValue>
	class Score<TScoreValue, ExtendedScore>
	{
		public:
		TScoreValue matchScore_;
		TScoreValue mismatchScore_;
		TScoreValue gapOpenInsertScore_;
		TScoreValue gapExtendInsertScore_;
		TScoreValue gapOpenDeleteScore_;
		TScoreValue gapExtendDeleteScore_;

		Score() :
		matchScore_(0), mismatchScore_(0), gapOpenInsertScore_(0), gapExtendInsertScore_(0), gapOpenDeleteScore_(0), gapExtendDeleteScore_(0)
		{
			SEQAN_CHECKPOINT;
		}

		Score(TScoreValue m, TScoreValue mm, TScoreValue goi, TScoreValue gei, TScoreValue god, TScoreValue ged) :
		 matchScore_(m), mismatchScore_(mm), gapOpenInsertScore_(goi), gapExtendInsertScore_(gei), gapOpenDeleteScore_(god), gapExtendDeleteScore_(ged)
		{
			SEQAN_CHECKPOINT;
		}

		Score(Score const & other) :
			matchScore_(other.matchScore_), mismatchScore_(other.mismatchScore_), gapOpenInsertScore_(other.gapOpenInsertScore_), gapExtendInsertScore_(other.gapExtendInsertScore_), gapOpenDeleteScore_(other.gapOpenDeleteScore_), gapExtendDeleteScore_(other.gapExtendDeleteScore_)
		{
			SEQAN_CHECKPOINT;
		}
		 ~Score()
		 {
			 SEQAN_CHECKPOINT;
		 }
	};

	template <typename TValue>
	inline TValue
		scoreMatch(Score<TValue, ExtendedScore> const & me)
	{
		SEQAN_CHECKPOINT;
		return me.matchScore_;
	}

	template <typename TValue>
	inline TValue
		scoreMismatch(Score<TValue, ExtendedScore> const & me)
	{
		SEQAN_CHECKPOINT;
		return me.mismatchScore_;
	}

	template <typename TValue>
	inline TValue
		scoreGapOpenInsert(Score<TValue, ExtendedScore> const & me)
	{
		SEQAN_CHECKPOINT;
		return me.gapOpenInsertScore_;
	}

	template <typename TValue>
	inline TValue
		scoreGapExtendInsert(Score<TValue, ExtendedScore> const & me)
	{
		SEQAN_CHECKPOINT;
		return me.gapExtendInsertScore_;
	}

	template <typename TValue>
	inline TValue
		scoreGapOpenDelete(Score<TValue, ExtendedScore> const & me)
	{
		SEQAN_CHECKPOINT;
		return me.gapOpenDeleteScore_;
	}

	template <typename TValue>
	inline TValue
		scoreGapExtendDelete(Score<TValue, ExtendedScore> const & me)
	{
		SEQAN_CHECKPOINT;
		return me.gapExtendDeleteScore_;
	}

	template <typename TValue, typename TPos1, typename TPos2, typename TStr1, typename TStr2>
	inline TValue
	scoreGapOpenVertical(Score<TValue, ExtendedScore> const & me,
						 TPos1,
						 TPos2,
						 TStr1 const &,
						 TStr2 const &)
	 {
		SEQAN_CHECKPOINT;
		return scoreGapOpenDelete(me);
	 }

	template <typename TValue, typename TPos1, typename TPos2, typename TStr1, typename TStr2>
	inline TValue
	scoreGapExtendVertical(Score<TValue, ExtendedScore> const & me,
						 TPos1,
						 TPos2,
						 TStr1 const &,
						 TStr2 const &)
	 {
		SEQAN_CHECKPOINT;
		return scoreGapExtendDelete(me);
	 }

	template <typename TValue, typename TPos1, typename TPos2, typename TStr1, typename TStr2>
	inline TValue
	scoreGapOpenHorizontal(Score<TValue, ExtendedScore> const & me,
						 TPos1,
						 TPos2,
						 TStr1 const &,
						 TStr2 const &)
	 {
		SEQAN_CHECKPOINT;
		return scoreGapOpenInsert(me);
	 }

	template <typename TValue, typename TPos1, typename TPos2, typename TStr1, typename TStr2>
	inline TValue
	scoreGapExtendHorizontal(Score<TValue, ExtendedScore> const & me,
						 TPos1,
						 TPos2,
						 TStr1 const &,
						 TStr2 const &)
	 {
		SEQAN_CHECKPOINT;
		return scoreGapExtendInsert(me);
	 }
}

#endif /* SCORE_JOURNAL_H_ */
