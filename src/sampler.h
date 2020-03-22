#ifndef SAMPLER_H
#define SAMPLER_H
#include "types.h"

inline std::vector<Vec2> generateJitteredSamplesOnUnitSquare(RandomCtx& rctx, uint32_t totalSamples)
{
	const uint32_t sampleCount = sqrt(totalSamples);
	std::vector<Vec2> out = {};
	out.reserve(sampleCount * sampleCount);
	
	for(uint32_t i = 0; i < sampleCount; i++)
	{
		for(uint32_t j = 0; j < sampleCount; j++)
		{
			Vec2 sample = {};
			sample.u = (j + rctx.distr(rctx.rengine)) / sampleCount; 
			sample.v = (i + rctx.distr(rctx.rengine)) / sampleCount;
			out.push_back(sample);
		}
	}

	return out;
}

inline std::vector<Vec2> mapJitteredSquareSamplesToUnitCircle(const std::vector<Vec2>& squareSamples)
{
	std::vector<Vec2> out = {};
	out.resize(squareSamples.size());

	for(uint32_t i = 0; i < squareSamples.size(); i++)
	{
		float r = 0;
		float phi = 0;

		//map square sample from [0,1] to [-1, 1] 2d range
		float x = squareSamples[i].u * 2.f - 1.f;
		float y = squareSamples[i].v * 2.f - 1.f;

		if(x > -y) // 1, 2 sectors
		{
			if(x > y)
			{
				r = x;
				phi = y / x;
			}
			else
			{
				r = y;
				phi = 2 - x / y;
			}
		}
		else // 3,4 sectors
		{
			if(x < y)
			{
				r = -x;
				phi = 4 + y / x;
			}
			else
			{
				r = -y;
				if(std::abs(y) > EPSILON)
				{
					phi = 6 - x / y;
				}
				else
				{
					phi = 0;
				}
			}
		}
		
		phi *= PI / 4.f;
		out[i] = {r * cos(phi), r * sin(phi)};
	}

	return out;
}

inline Vec2 sampleUnitCircle(RandomCtx& rctx, const std::vector<Vec2>& unitSamples)
{
	const uint32_t randomIndex = rctx.distr(rctx.rengine) * (unitSamples.size() - 1);
	return unitSamples[randomIndex];
}

#endif